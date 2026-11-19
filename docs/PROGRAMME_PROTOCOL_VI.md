# Protocol chương trình sắn

Phục vụ bài review sắn riêng. Public protocols và primary studies đã đọc được ánh xạ trong `references/cassava_sources.json`; giới hạn đọc và truy cập giữ nguyên.

## Sản phẩm và di truyền

Sản phẩm là clone có năng suất củ tươi, dry matter/starch, sức khỏe và chất lượng phù hợp mục đích sử dụng. Parent recycling hướng tới phần truyền được; clonal advancement khai thác total genetic value. Một clone xếp hạng cao về A+D không tự là bố mẹ có breeding value tốt nhất.

AM560-2 v8 có accession **GCF_001659605.2/GCA_001659605.2**; v6 dùng accession kết thúc `.1` và được giữ như phiên bản lịch sử. TME204 H1 **GCA_020916445.1** và H2 **GCA_020916425.1** là hai haplotype của một cultivar lưỡng bội. Không gọi chúng là hai giống hoặc một cây tứ bội. Cây sắn có 18 chromosome hạt nhân trong bộ đơn bội; số đếm có thêm plastid/mitochondrion không được đổi thành 20 chromosome hạt nhân.

Pilot dùng 18 chromosome tổng hợp, phased LD blocks, QTLAD và marker tách vị trí. Morgan/map lengths, QTL effects, dominance, allele frequencies và founder structure là giả định. Không lấy assembly sequence làm dữ liệu genotype Việt Nam hoặc bản đồ tái tổ hợp. Biến `gv1` là total genetic value trên thang cohort ban đầu, không phải additive breeding value.

## Thiết kế chọn giống

1. Khóa product profile, clone/parent ID, pedigree và fingerprints. Ghi riêng nguồn botanical seed và vật liệu vô tính.
2. Chọn bố mẹ từ mô hình additive riêng, kiểm tra relatedness/flowering/crossing capacity; tạo progeny bằng meiosis. Code hiện dùng family cap, chưa tối ưu contribution hoặc flowering success.
3. Nhân clone và khảo nghiệm theo các giai đoạn: vườn cây con → CET/SRT → PYT → AYT → UYT/đa môi trường → đánh giá nông hộ/chế biến. Pilot gộp quan sát vào ba địa điểm, ba chỉ số năm/vụ, hai chế độ quản lý và hai block lặp; nó chưa triển khai riêng từng giai đoạn hoặc giới hạn nguồn hom.
4. Đo củ tươi, dry matter, starch/quality và disease/pathogen bằng phương pháp tham chiếu phù hợp. Parent selection và clone selection được ghi thành quyết định khác nhau.
5. Chuyển dòng qua đánh giá và hệ thống giống có identity/health/traceability. Hồ sơ sản phẩm, sức khỏe vật liệu và đường công bố/lưu hành tại Việt Nam được đối chiếu trong [workflow vận hành](CASSAVA_WORKFLOW_VI.md); không dùng ngưỡng mô phỏng làm tiêu chuẩn chứng nhận.

Chu kỳ 5 năm trong config là lịch giả định. Biện pháp nhân nhanh, tỷ lệ sống, seed set và mức rút ngắn chu kỳ chưa được hiệu chuẩn. Nền tảng IITA/CIAT/NextGen và BASICS phục vụ thiết kế nhánh sắn; SOP ngô lai của doanh nghiệp không được gán cho chương trình sắn.

## UAV và kiểu hình

UAV cung cấp tán, hình thái, sắc tố và nhiệt tại các thời điểm đã khóa trước quyết định. Nó không chụp trực tiếp củ dưới đất. Chụp ảnh củ đã thu hoạch, spectroscopy trên củ nguyên/nghiền, radar và aerial canopy imaging có các target/reference khác nhau.

Nghiên cứu UAV Selvaraj2020 có held-out field replication; đó không phải một kiểm định trên genotype/year/site mới. Nghiên cứu ảnh hình dạng củ Yonis2020 dùng củ đã thu hoạch. Cả hai không xác nhận pipeline hiện tại ngoài đồng; chúng xác định những phép so sánh đáng thử và giới hạn diễn giải.

Code hiện tạo early canopy features có nhiễu và missingness. Chuỗi thực gồm raw flights, calibration/geometric QC, plot interior/segmentation, time series và harvest/lab references. Whole-plot, clone, family, site, season, flight và image-source overlap phải được kiểm soát khi chia train/test.

## Dự đoán và quyết định

So sánh số đo đồng ruộng cũ, genomic và genomic+UAV bằng cùng số vật liệu/field records. Tần số allele, imputation và centering chỉ fit training. Kiểm định harvest mới, actual parents mới, site mới hoặc đồng thời cả ba. Quan hệ họ hàng xa có thể còn ngay cả khi bố mẹ không trùng.

Xếp hạng clone dùng ba predicted traits và trọng số giả định. Parent index dùng additive-only model riêng, không lấy conditional marker coefficient trong mô hình UAV làm toàn bộ breeding value. Vì hai nhánh genomic dùng cùng parent predictor, pilot chỉ kiểm tra giá trị UAV bổ sung đối với forecast/ranking sản phẩm. Nó không chứng minh thêm parent genetic response hay giải bài toán truyền bệnh qua hom.

## Nhân giống sạch và quản lý

Quy trình đề xuất theo dõi mother stock, identity, pathogen testing, propagation batch, health inspections, loại bỏ khi có căn cứ, phân phối và kết quả ngoài ruộng. Ảnh chẩn đoán phải được đối chiếu với tác nhân và có thể bỏ sót cây không triệu chứng. Clean-seed multiplication và epidemic được giữ như nhánh dữ liệu/mô hình riêng; pilot không mô phỏng chúng.

Thí nghiệm policy hiện phân treatment ngẫu nhiên, có test fields độc lập và hai potential outcomes tổng hợp. Sensor/proxy có nhiễu, cost chuyển thành harvest equivalent là giả định. Kết quả không phải liều phân/tưới/thuốc, mức nước tiết kiệm, prevalence hoặc hiệu quả kiểm soát SLCMV đã đo.

## Căn cứ và hiệu chuẩn tiếp

Mười nguồn ở ledger gồm registry/genome, genetic prediction, root/UAV phenotyping và các hệ thống clean stock. Không coi tất cả là full-paper/supplement assessment. Các thông số pilot không được suy ra từ một company announcement hoặc từ accuracy trong một thí nghiệm khác.

Cần hiệu chuẩn founder/maps, LD, genetic/residual covariance, target/reference measurements, genotype-by-environment, sensor error, disease pressure, propagation/crossing success, capacity, economic weights và latency trước khi suy cho chương trình Việt Nam. Mô phỏng hiện tại là thí nghiệm thiết kế đồng hành review; synthetic results không đi vào inventory272 nghiên cứu hoặc các phép so sánh thực nghiệm của bản thảo.
