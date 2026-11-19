# Phương trình di truyền

Đối chiếu mã ngày 19/11/2026: [genetics.R](../src/uavbreed/genetics.R), [models.py](../src/uavbreed/models.py), [programme.py](../src/uavbreed/programme.py). Các phương trình dưới đây mô tả nhánh sắn. Các hệ số sinh dữ liệu là giả định mô phỏng.

## 1. Haplotype, liên kết và giảm phân

Hai haplotype lưu pha liên kết; giao tử nhận các đoạn qua giảm phân và tái tổ hợp. Mã khởi tạo khối allele với xác suất làm mới \(1-e^{-\Delta m/\ell}\), trong đó \(m,\ell\) có đơn vị Morgan. Đây là LD tổng hợp, không phải lịch sử quần thể đã ước lượng. Bản đồ có khoảng cách Morgan; assembly cung cấp tọa độ bp. Không suy bản đồ tái tổ hợp từ chiều dài assembly. [AlphaSimR: newMapPop, Details](https://gaynorr.github.io/AlphaSimR/reference/newMapPop.html).


## 2. Giá trị di truyền và quan sát ngoài đồng

Cấu trúc tổng quát là \(Y=\mu+G+E+G\times E+f(M)+B+\epsilon\). Trong `field_trial()`, \(g_k=\texttt{gvk}\) đã chuẩn hóa theo quần thể chu kỳ 0. Biến mã `a` chính là \(g_1\), không phải breeding value additive. Công thức năng suất đang dùng:

\[
\begin{aligned}
d&=\operatorname{logit}^{-1}(2p-0.65g_4-1),\\
\eta&=\mu_y+s_y[g_1+0.30g_3w+0.20g_1n-0.80w+0.30n-0.50d],\\
Y&=\eta+s_y(B+\epsilon).
\end{aligned}
\]

\(w,n,p\) là `water`, `n_supply`, `pressure` tổng hợp; quản lý thay đổi \(w,n\). Không có hiệu ứng quản lý ngẫu nhiên riêng. Chất lượng dùng \(\mu_q+s_q(g_2+0.10w+\epsilon_q)\); bệnh dùng \(\operatorname{clip}(d+\epsilon_d,0,1)\). UAV là proxy nhiễu của di truyền/môi trường, chưa phải ảnh hoặc quá trình bay vật lý. Phân biệt cơ chế sinh dữ liệu với mô hình ước lượng theo [Morris et al., §3.2–3.3](https://doi.org/10.1002/sim.8086).

## 3. Hiệu ứng chức năng và breeding value

Mô hình AD ở locus có giá trị chức năng \((-a,d,a)\) cho dosage \((0,1,2)\). Dưới HWE, với tần số allele được đếm là \(p\), \(q=1-p\):

\[
\alpha=a+(q-p)d,\quad A_i=\sum_l(x_{il}-2p_l)\alpha_l,\quad
G_i=\mu_G+A_i+D_i.
\]

Dominance chức năng góp vào \(\alpha\); additive không đồng nghĩa chỉ cộng các \(a_l\). Đẳng thức \(V_G=V_A+V_D\) cần điều kiện trực giao; LD và lệch HWE có thể thêm covariance. [Vitezica et al., Theory, (2)–(4)](https://doi.org/10.1534/genetics.113.155176). Mã xuất `gv*` và `bv*` riêng; `genetic_sd` hiện là SD của `gv1`, không phải \(\sigma_A\). [AlphaSimR: bv](https://gaynorr.github.io/AlphaSimR/reference/bv.html).

## 4. Kernel diploid của clone

Nhánh không-ngô dùng \(Z=x-2p\), \(W_{x=0,1,2}=(-2p^2,2pq,-2q^2)\):

\[
K_A=ZZ'/\sum_l2p_lq_l,\qquad K_D=WW'/\sum_l(2p_lq_l)^2.
\]

Đây là mã hóa HWE [Vitezica et al., Theory](https://doi.org/10.1534/genetics.113.155176). Clone chọn lọc, quan hệ gia đình và dị hợp bất thường không bảo đảm trực giao A/D. Với sắn, giữ nguyên clone khai thác \(G\); chọn bố mẹ hướng tới phần truyền được. Hai mục tiêu khác nhau, như [Wolfe et al., Discussion](https://doi.org/10.1534/g3.116.033332); nguồn này không hiệu chuẩn phương sai cho repo.

## 5. Ridge và hai mục tiêu dự báo

`fit_ridge` giải:

\[
\min_{B,b_0}\|Y-\mathbf1b_0'-XB\|_F^2+\lambda\|B\|_F^2,\qquad\lambda=1.
\]

Đây là [Ridge, objective/alpha](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html), không phải REML. Median, chỉ báo thiếu và centering đều học từ training; `with_std=False` không chia SD. Cùng penalty áp dụng lên các cột sau xử lý, nên thang đo ảnh hưởng shrinkage. Ba đầu ra được hồi quy song song; không ước lượng phương sai A/D hoặc covariance đa tính trạng.

Sản phẩm dùng \(X=[\Phi_A,\Phi_D,E]\), thêm UAV ở `genomic_uav`. Mô hình bố mẹ được fit **riêng**, \(X_P=[\Phi_A,E]\), bỏ dominance và UAV; điểm marker là \(\widehat A^{marker}=\Phi_A\widehat B_A\). Đây vẫn là thành phần marker fit từ quan sát, không phải `bv` thật. Hai chiến lược genomic dùng cùng mô hình bố mẹ; lợi ích UAV hiện nằm ở xếp hạng sản phẩm.

## 6. Chỉ số nhiều tính trạng

Với năng suất, chất lượng, bệnh:

\[
I=\sum_{j=1}^3w_j\frac{\widehat Y_j-\bar Y_{j,tr}}{s_{j,tr}},\quad
w=(0.65,0.20,-0.15);\qquad
I_P=\sum_jw_j\widehat A^{marker}_j/s_{j,tr}.
\]

Trọng số là lựa chọn kịch bản; chưa phải chỉ số Smith–Hazel tối ưu từ covariance và giá trị kinh tế. [Hazel, Introduction và xây dựng chỉ số](https://doi.org/10.1093/genetics/28.6.476). Mã lấy trung bình điểm theo genotype. Bố mẹ sắn có trần hai cá thể trên mỗi nhãn gia đình; nhãn gia đình chưa đo mức kinship giữa các bố mẹ.

## 7. Gain theo năm và bất định

Chuẩn lý thuyết là:

\[
\Delta A/\text{năm}\approx i\,r_{I,A}\sigma_A/L,\qquad
Q=\frac{i_Ur_U\sigma_{A,U}/L_U}{i_Br_B\sigma_{A,B}/L_B}.
\]

\(i\) là chênh lệch điểm chọn lọc chuẩn hóa, \(\sigma_A=\operatorname{SD}(A)\), \(L\) là số năm/vòng tái tổ hợp–đánh giá–chọn bố mẹ; U/B chỉ chiến lược UAV/nền. Cần xác định quần thể và mục tiêu, hồi quy \(A\mid I\) gần tuyến tính; các đại lượng không tự giữ cố định qua nhiều vòng. [Rutkoski, §4](https://doi.org/10.1016/bs.agron.2019.05.001), [Bảng 2 và §4.5 trong bản tác giả](https://smallgrains.web.illinois.edu/wp/wp-content/uploads/2022/09/A-practical-guide-to-genetic-gain-my-copy-2.pdf). \(r_{I,A}\) là Pearson với breeding value, không phải Spearman với `yield_obs`. Ràng buộc gia đình cần được đánh giá bằng mô phỏng; không mặc định công thức một quần thể áp dụng nguyên vẹn.

Mã chưa tính công thức gain này. `selected_potential_yield` là kết quả sản phẩm; `mean_genetic_value` là trung bình \(G\). Calendar \(cL\) giả định cùng \(L\) giữa chiến lược, với vật liệu đánh giá đã có ở chu kỳ 0; không đo rút ngắn chu kỳ hay thời gian phát hành. Tỷ số \(Q\) không ổn định khi mẫu số gần 0.

Với chênh lệch ghép cặp mỗi replicate \(d_r\):

\[
\bar d=R^{-1}\sum_rd_r,\qquad \mathrm{MCSE}(\bar d)=s_d/\sqrt R.
\]

Đơn vị lặp độc lập là replicate, không phải ô/thửa. Phân vị giữa replicate khác khoảng tin cậy của trung bình. Pilot có 8 replicate; công thức MCSE và khoảng t của trung bình được tính từ các CSV trong `report.py`; `programme.py` xuất từng lần lặp. MCSE không bao gồm bất định do tham số ngoài đồng chưa hiệu chuẩn. [Morris et al., §4.3, §5](https://doi.org/10.1002/sim.8086). Thí nghiệm quản lý ngẫu nhiên đánh giá lợi ích can thiệp sau chi phí, tách khỏi gain di truyền.
