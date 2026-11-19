# Triển khai giống đã chọn vào thửa quản lý mô phỏng

Chạy riêng, không gọi lại engine di truyền:

```powershell
analysis/env/python/Scripts/python.exe -I -X utf8 run_deployment.py
```

Mặc định đọc `analysis/results/validated_pilot`, ghi thư mục mới `analysis/results/genotype_deployment`. `--source` và `--output` cho phép chỉ định đường dẫn khác. Thư mục đầu ra phải trống và nằm ngoài nguồn đóng băng.

## Dòng dữ liệu

Điểm `genomic_uav` của chu kỳ cuối được lấy trung bình theo ID. Cà phê giữ 6 clone đầu bảng và kiểm tra toàn bộ ma trận tương hợp đã lưu với allele S của chúng. Ngô giữ `floor(0.20 × số tổ hợp)` tổ hợp đầu bảng; nhánh `cassava` áp dụng cùng tỷ lệ cho clone khi dùng với repo sắn. Chọn bằng điểm sản phẩm, không dùng `gv` để chọn.

Mỗi replicate tạo riêng 600 thửa huấn luyện và 1.000 thửa kiểm tra theo cấu hình nguồn. ID được phân bổ cân bằng trong từng tập, nên đây là các thửa độc lập của **những giống đã biết**, không phải thử nghiệm giống mới. Các chính sách dùng chung thửa kiểm tra. Ngẫu nhiên hóa xử lý được phân tầng theo ID trong từng tập.

## Phương trình đang chạy

Với thửa \(f\) trồng giống \(g\), \(s_f\sim U(0,1)\), \(v_g=\texttt{gv1}\), \(r_g=\texttt{gv3}\). Giữ đơn vị năng suất và thang \(k\) của cấu hình cây trồng:

\[
Y_f(0)=\mu+k(v_g+0.20r_g-s_f),\qquad
\delta_f=0.90s_f+0.15r_g-0.15.
\]

Với xác suất `nonresponse_rate`, thay \(\delta_f\) bằng \(-0.12\). Sau đó \(Y_f(1)=Y_f(0)+k\delta_f\). Các hệ số này tái sử dụng ý tưởng DGP quản lý tổng hợp, không phải ước lượng từ nông hộ hay công ty.

\[
S_f^{obs}=s_f+N(0,\sigma_{sensor}^2),\quad
R_f^{obs}=r_g+N(0,0.40^2),\quad
Y_f^{obs}=Y_f(T_f)+N(0,(0.35k)^2).
\]

Đầu vào mô hình chỉ gồm one-hot ID giống, hai proxy quan sát và một biến nhiễu quan sát độc lập. Hai hồi quy `HistGradientBoostingRegressor` được fit trên thửa huấn luyện theo \(T=0/1\). Với chi phí \(c=k\times\texttt{treatment_cost_yield_units}\), xử lý khi \(\widehat Y(1)-\widehat Y(0)>c\). Hệ số, kích thước dữ liệu, seed và thiết lập learner được lưu trong `deployment_config.json`.

So sánh `none`, `all`, quy tắc sensor > 0.5, chính sách học và oracle trên cùng thửa. Tiện ích là \(Y_f(a_f)-ca_f\). Oracle dùng hai kết quả tiềm năng thật **sau khi fit**, chỉ làm cận trên. Mô hình không nhận `gv1`, `gv3`, stress thật, trạng thái không đáp ứng hoặc kết quả tiềm năng. Việc tách DGP, dữ liệu quan sát và mục tiêu đánh giá theo [Morris et al., §3.2–3.3 và §4.3](https://doi.org/10.1002/sim.8086).

## Kết quả và phạm vi

- `policy_metrics.csv`: kết quả theo replicate/cây/chính sách, tỷ lệ xử lý, lợi ích sau chi phí và chênh lệch với không xử lý.
- `selected_candidates_TRUTH_ONLY.csv`: ID được chọn, thứ hạng và giá trị di truyền dùng sinh dữ liệu.
- `observations_SYNTHETIC.csv`: dữ liệu quan sát, ID và tập train/test.
- `field_truth_TRUTH_ONLY.csv`, `test_field_policies_TRUTH_ONLY.csv`: ID thửa/giống, hai kết quả tiềm năng, quyết định và tiện ích của từng chính sách.
- `coffee_selected_compatibility.csv`: ma trận S của bộ 6 clone; không đưa tương hợp vào hàm năng suất, không mô phỏng phối trồng hoặc đậu quả.
- Các tệp config/manifest ghi phiên bản, seed, hash đầu vào, mã triển khai và đầu ra. Nguồn sinh học trong `references/*_sources.json` hỗ trợ bối cảnh; không hiệu chuẩn DGP này.

Đây là triển khai giả định của giống đã chọn vào thửa mới cùng phân bố tổng hợp. Kết quả không đo khả năng chuyển sang giống mới, môi trường sản xuất khác, đậu quả cà phê, hay phản hồi từ nông hộ. Pipeline chọn giống và các kết quả đã đánh giá trước đó được giữ nguyên.
