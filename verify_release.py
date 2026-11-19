"""Verify release files and the numerical sources of Table S6."""
import hashlib
import json
from pathlib import Path
import tempfile
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((ROOT / 'SHA256.json').read_text(encoding='utf-8'))
    for name, expected in manifest.items():
        path = (ROOT / name).resolve()
        assert path.is_relative_to(ROOT) and path.is_file(), name
        assert sha(path) == expected, name
    pilot = ROOT / 'analysis/results/validated_pilot'
    run = json.loads((pilot / 'run_manifest.json').read_text(encoding='utf-8'))
    for key, name in [('config_sha256', 'configs/pilot.json'),
                      ('engine_sha256', 'src/uavbreed/genetics.R'),
                      ('programme_sha256', 'src/uavbreed/programme.py')]:
        assert run[key] == sha(ROOT / name), key
    assert run['replicates'] == 8 and run['cycles'] == 2 and run['crops'] == ['cassava']
    deployment = ROOT / 'analysis/results/genotype_deployment'
    provenance = json.loads((deployment / 'source_manifest.json').read_text(encoding='utf-8'))
    for record in provenance:
        path = (ROOT / record['path']).resolve()
        assert path.is_relative_to(ROOT) and sha(path) == record['sha256'], record['path']
    for name, expected in json.loads((deployment / 'run_manifest.json').read_text(encoding='utf-8'))['outputs'].items():
        assert sha(deployment / name) == expected, name
    with tempfile.TemporaryDirectory(prefix='cassava_s6_') as temporary:
        output = Path(temporary) / 's6'
        subprocess.run([sys.executable, '-I', '-X', 'utf8', str(ROOT / 'scripts/reproduce_s6.py'),
                        '--programme-root', str(ROOT), '--output', str(output)], check=True, cwd=ROOT)
        for name in ['Supplementary_Table_S6_Simulation_Contrasts.csv',
                     'Supplementary_Table_S6_Simulation_Replicates.csv']:
            assert (output / name).read_bytes() == (ROOT / 'data/expected_s6' / name).read_bytes(), name
    print(json.dumps({'release_files_checked': len(manifest), 'source_hashes_match': True,
                      'seven_contrasts_and_56_values_match': True}))


if __name__ == '__main__':
    main()
