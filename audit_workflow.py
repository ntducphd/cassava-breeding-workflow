"""Check the cassava pilot and export the operational-workflow assessment."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
COVERAGE = [
    ('W01', 'Product profile and target environments', 'partial', 'Market, target region, quality and cost thresholds', 'SCHEMES2021'),
    ('W02', 'Germplasm identity and entry health', 'not_modelled', 'Fingerprint, accession, pathogen panel and material provenance', 'SEED_TOOLKIT'),
    ('W03', 'Crossing and botanical seed', 'partial', 'Meiosis is simulated; flowering, seed set and germination require records', 'CIAT2016'),
    ('W04', 'Seedling nursery', 'not_modelled', 'Survival, seedling identity and available clonal planting material', 'CIAT2016'),
    ('W05', 'CET or single-row trial', 'not_separate', 'Clone plots, check varieties, harvest age and stem availability', 'SCHEMES2021;IITA06'),
    ('W06', 'Preliminary yield trial', 'not_separate', 'Replicated trial and product-specific advancement criteria', 'SCHEMES2021;IITA06'),
    ('W07', 'Advanced yield trial', 'not_separate', 'Connected sites/seasons, trial quality and GxE analysis', 'SCHEMES2021;IITA06'),
    ('W08', 'Uniform or regional yield trial', 'not_separate', 'Target-environment performance, checks and uncertainty', 'SCHEMES2021;IITA06'),
    ('W09', 'On-farm and processor evaluation', 'not_modelled', 'Observed farm outcomes, usability and processor acceptance', 'CIAT2016'),
    ('W10', 'Parent versus clonal selection', 'implemented_synthetic', 'Additive parent target versus clonal product performance', 'CIAT2016'),
    ('W11', 'UAV and root reference measurements', 'partial', 'Directly generated proxies require flight and harvest calibration', 'IITA06'),
    ('W12', 'Stem health and multiplication', 'not_modelled', 'Mother stock, laboratory results, propagation batches and lot tracing', 'SEED_TOOLKIT;CMD2_TC2018;CMD2_2022'),
    ('W13', 'Vietnam information disclosure', 'external_requirement', 'Applicable TCVN/TCCS, variety information and cultivation process', 'VN_TT12;VN_NQ66;VN_TT15'),
    ('W14', 'DUS and plant variety rights', 'external_requirement', 'Separate route where protection or examination applies', 'UPOV321'),
    ('W15', 'Budget, calendar and sensitivity', 'not_calibrated', 'Stage capacity, actual elapsed time, costs and sensitivity analysis', 'SCHEMES2021'),
    ('W16', 'Selected-clone management', 'implemented_synthetic', 'Independent fields of the same known selected clone IDs', ''),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'analysis/validation/workflow_20261119')
    args = parser.parse_args()
    output = args.output.resolve()
    cfg = json.loads((ROOT / 'configs/pilot.json').read_text(encoding='utf-8'))
    saved = ROOT / 'analysis/results/validated_pilot'
    assert cfg == json.loads((saved / 'config_used.json').read_text(encoding='utf-8'))
    assert set(cfg['crops']) == {'cassava'}
    assert cfg['crops']['cassava']['chromosomes'] == 18
    splits = pd.read_csv(saved / 'split_audit.csv')
    assert len(splits) == cfg['replicates'] * len(cfg['strategies']) * cfg['cycles'] * 4
    assert splits.row_overlap.eq(0).all()
    assert splits.train_year_max.lt(splits.test_year_min).all()
    held = splits.validation.isin(['new_parents', 'new_parents_site'])
    assert splits.loc[held, 'all_parent_overlap'].eq(0).all()
    assert splits.loc[held, 'genotype_overlap'].eq(0).all()
    audits = list(saved.glob('snapshots/**/genetics_audit.json'))
    for path in audits:
        audit = json.loads(path.read_text(encoding='utf-8'))
        assert audit['ploidy'] == 2 and audit['chromosomes'] == 18 and audit['qtl_marker_overlap'] == 0
    assert audits
    selected = pd.read_csv(ROOT / 'analysis/results/genotype_deployment/selected_candidates_TRUTH_ONLY.csv',
                           dtype={'genotype_id': str})
    expected = max(1, int(cfg['genetics']['families'] * cfg['genetics']['progeny_per_family'] * cfg['selection_fraction']))
    counts = selected.groupby(['replicate', 'crop']).genotype_id.nunique()
    assert len(counts) == cfg['replicates'] and counts.eq(expected).all()
    ledger_path = ROOT / 'references/workflow_sources_20261119.json'
    ledger = json.loads(ledger_path.read_text(encoding='utf-8'))
    sources = {item['id']: item for item in ledger}
    for row in COVERAGE:
        assert all(ident in sources for ident in row[4].split(';') if ident)
    external_snapshots = []
    for item in ledger:
        if item.get('downloaded') and item.get('usable', True):
            snapshot = ROOT / item['file']
            if snapshot.is_file():
                assert sha(snapshot) == item['sha256'], item['id']
            else:
                external_snapshots.append(item['id'])
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'workflow_coverage.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['ID', 'Stage', 'PilotCoverage', 'RequiredForOperationalUse', 'EvidenceIDs'])
        writer.writerows(COVERAGE)
    source_code = list((ROOT / 'src/uavbreed').glob('*.py')) + [ROOT / 'src/uavbreed/genetics.R', ROOT / 'configs/pilot.json',
                  ROOT / 'docs/CASSAVA_WORKFLOW_VI.md', ROOT / 'docs/GENETIC_EQUATIONS.md',
                  ROOT / 'docs/PROGRAMME_PROTOCOL_VI.md', ledger_path,
                  ROOT / 'tests/test_deployment.py', ROOT / 'tests/test_scientific_integrity.py']
    report = {'assessment_date': '2026-11-19', 'simulation_integrity_passed': True,
              'founder_replicates': cfg['replicates'], 'evaluated_cycles': cfg['cycles'],
              'genetics_audits_checked': len(audits), 'split_rows_checked': len(splits),
              'parent_and_genotype_overlap_in_new_parent_tests': 0, 'selected_clones_per_replicate': expected,
              'operational_stages_assessed': len(COVERAGE), 'field_programme_validated': False,
              'empirical_calibration_complete': False, 'variety_release_authorized': False,
              'unused_timing_parameter': 'phenotyping.early_fraction_of_season',
              'manual_stage_assessment': 'docs/CASSAVA_WORKFLOW_VI.md',
              'source_sha256': {p.relative_to(ROOT).as_posix(): sha(p) for p in source_code},
              'reference_records': len(ledger), 'external_source_snapshots_not_bundled': external_snapshots}
    (output / 'workflow_checks.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
