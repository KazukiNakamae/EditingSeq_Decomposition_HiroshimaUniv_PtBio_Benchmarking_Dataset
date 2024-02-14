import numpy as np
import pandas as pd
import glob
import os
import json
import openpyxl
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats
import seaborn as sns
import datetime

from analysis_library import SynDNAData
from analysis_library import analyze_one_condition_data
from analysis_library import analyze_one_condition_data_custom
from analysis_library import comp_accuracies
from analysis_library import comp_two_factor

dt_stamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
res_dir = 'result_' + dt_stamp
os.makedirs(res_dir)

### 個別解析
fwd_samples = [str(i) for i in list(range(1, 9, 1))] + ['18'] + ['9']
rev_samples = [str(i) for i in list(range(11, 16, 1))] + ['19'] + ['16']
low_samples = ['5', '22', '23']

# backstep/L1正則化なし
backstep_nonl1_patt = 'PtWAVE_*_on_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(backstep_nonl1_patt):
  print(name)
# 問題なし
backstep_nonl1_data = SynDNAData(patt=backstep_nonl1_patt, max_deletion=100)
analyze_one_condition_data(backstep_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_backstep_nonl1_PtWAVE")

# backstep/L1正則化あり
backstep_l1_patt = 'PtWAVE_*_on_L1_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(backstep_l1_patt):
  print(name)
# 問題なし
backstep_l1_data = SynDNAData(patt=backstep_l1_patt, max_deletion=100)
analyze_one_condition_data(backstep_l1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_backstep_l1_PtWAVE")

# all/L1正則化なし
all_nonl1_patt = 'PtWAVE_*_on_all_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(all_nonl1_patt):
  print(name)
# 問題なし
all_nonl1_data = SynDNAData(patt=all_nonl1_patt, max_deletion=100)
analyze_one_condition_data(all_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_all_nonl1_PtWAVE")

# all/L1正則化あり
all_l1_patt = 'PtWAVE_*_on_L1_all_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(all_l1_patt):
  print(name)
# 問題なし
all_l1_data = SynDNAData(patt=all_l1_patt, max_deletion=100)
analyze_one_condition_data(all_l1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_all_l1_PtWAVE")

# random/L1正則化なし
random_nonl1_patt = 'PtWAVE_*_on_random_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(random_nonl1_patt):
  print(name)
# 問題なし
random_nonl1_data = SynDNAData(patt=random_nonl1_patt, max_deletion=100)
analyze_one_condition_data(random_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_random_nonl1_PtWAVE")

# random/L1正則化あり
random_l1_patt = 'PtWAVE_*_on_L1_random_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(random_l1_patt):
  print(name)
# 問題なし
random_l1_data = SynDNAData(patt=random_l1_patt, max_deletion=100)
analyze_one_condition_data(random_l1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_random_l1_PtWAVE")

# all/L1正則化なし/150bp欠失まで
seventyfive_all_nonl1_patt = 'PtWAVE_*_on_seventyfive_all_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(seventyfive_all_nonl1_patt):
  print(name)
# 問題なし
seventyfive_all_nonl1_data = SynDNAData(patt=seventyfive_all_nonl1_patt, max_deletion=150)
analyze_one_condition_data(seventyfive_all_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_seventyfive_all_nonl1_PtWAVE")

# backstep/L1正則化なし/150bp欠失まで
seventyfive_backstep_nonl1_patt = 'PtWAVE_*_on_seventyfive_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(seventyfive_backstep_nonl1_patt):
  print(name)
# 問題なし
seventyfive_backstep_nonl1_data = SynDNAData(patt=seventyfive_backstep_nonl1_patt, max_deletion=150)
analyze_one_condition_data(seventyfive_backstep_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_seventyfive_backstep_nonl1_PtWAVE")

# random/L1正則化なし/150bp欠失まで
seventyfive_random_nonl1_patt = 'PtWAVE_*_on_seventyfive_random_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(seventyfive_random_nonl1_patt):
  print(name)
# 問題なし
seventyfive_random_nonl1_data = SynDNAData(patt=seventyfive_random_nonl1_patt, max_deletion=150)
analyze_one_condition_data(seventyfive_random_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_seventyfive_random_nonl1_PtWAVE")

# all/L1正則化なし/20bp欠失まで
ten_all_nonl1_patt = 'PtWAVE_*_on_ten_all_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(ten_all_nonl1_patt):
  print(name)
# 問題なし
ten_all_nonl1_data = SynDNAData(patt=ten_all_nonl1_patt, max_deletion=100)
analyze_one_condition_data(ten_all_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_ten_all_nonl1_PtWAVE")

# backstep/L1正則化なし/20bp欠失まで
ten_backstep_nonl1_patt = 'PtWAVE_*_on_ten_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(ten_backstep_nonl1_patt):
  print(name)
# 問題なし
ten_backstep_nonl1_data = SynDNAData(patt=ten_backstep_nonl1_patt, max_deletion=100)
analyze_one_condition_data(ten_backstep_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_ten_backstep_nonl1_PtWAVE")

# random/L1正則化なし/20bp欠失まで
ten_random_nonl1_patt = 'PtWAVE_*_on_ten_random_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(ten_random_nonl1_patt):
  print(name)
# 問題なし
ten_random_nonl1_data = SynDNAData(patt=ten_random_nonl1_patt, max_deletion=100)
analyze_one_condition_data(ten_random_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_ten_random_nonl1_PtWAVE")
#!!!論文ポイント
# 定義するインデル範囲が狭すぎるとどんな手法でも精度が下がり、インデル効率推定の正解率も同程度の下がる
# したがってこの条件では既存のアルゴリズムと差別化できないので、インデル範囲が狭すぎた場合の結果について
# backstepのアルゴリズムの優位性を主張しないように気をつけなくてはならない

# all/L1正則化なし/150bp欠失/50bp挿入まで
seventyfiveINSfifty_all_nonl1_patt = 'PtWAVE_*_on_seventyfiveINSfifty_all_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(seventyfiveINSfifty_all_nonl1_patt):
  print(name)
# 問題なし
seventyfiveINSfifty_all_nonl1_data = SynDNAData(patt=seventyfiveINSfifty_all_nonl1_patt, max_deletion=150)
analyze_one_condition_data(seventyfiveINSfifty_all_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_seventyfiveINSfifty_all_nonl1_PtWAVE")
#!!!論文ポイント
# インデル種が一種しかない場合は、allアルゴリズムは極めて頑強であることが示された。
# したがってこの実験については、ラージデリーション検出がallとbackstepで可能であり
# BICの低下によるモデル安定化を図る上で、少なくともrandomのアルゴリズムは不適である
# という議論をすべきではないかと考えられる。
# ややポジショントークではあるが、この実験の意義は「変数選択による精度低下の懸念解消」
# という部分にフォーカスすべきかもしれない。


# backstep_nonl1_dataで、BIC値が二極化しているので分けて解析も行ってみる。
backstep_nonl1_data.contribs_df.query('BIC<10000')['-85']
backstep_nonl1_data.contribs_df.query('BIC<10000')['other_indels']
backstep_nonl1_data.contribs_df.query('BIC>10000')['-85']
backstep_nonl1_data.contribs_df.query('BIC>10000')['other_indels']
#!!!論文ポイント
# BICが低いものは-85を検出できておらず、インデル検出自体が低い傾向にある
# 考慮するインデルが減ればBICが小さくなるのは自明なので、当然の結果だが
# これはBICを意識しすぎると、インデルを低く見積もる可能性を示唆している
# この点はちゃんとディスカッションに書いておくべきである。


### PtWAVEパラメータ間解析
### 手法間比較
comp_accuracies(all_nonl1_data.contribs_df, random_nonl1_data.contribs_df, backstep_nonl1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_nonl1_100del3ins_comp')
#!!!論文ポイント
# backstepでは精度を維持しながらBICの低下をさせることに成功している
comp_accuracies(all_l1_data.contribs_df, random_l1_data.contribs_df, backstep_l1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_l1_100del3ins_comp')
comp_accuracies(seventyfive_all_nonl1_data.contribs_df, seventyfive_random_nonl1_data.contribs_df, seventyfive_backstep_nonl1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_seventyfive_nonl1_100del3ins_comp')
#!!!論文ポイント
# randomでは精度の低下が発生しているが、backstepでは精度を維持しながらBICの低下をさせることに成功している
comp_accuracies(ten_all_nonl1_data.contribs_df, ten_random_nonl1_data.contribs_df, ten_backstep_nonl1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_ten_nonl1_100del3ins_comp')
### L1 vs nonL1間比較
# backstep や allで見られる傾向として、L1はインデルが低濃度下ではインデル（-85またはediting eff）そのものを見逃す傾向があり、一方で高濃度下では正確にインデル（-85またはediting eff）を定量する傾向が見受けられる
# この点を期待インデル率20%以下の区画と期待インデル率80%以上の区画の散布図と決定係数を用いて明確に表現する
high_conc_samples_set=set(backstep_nonl1_data.contribs_df.query('expected_editing_eff>=80')['sampletype'].tolist())
low_conc_samples_set=set(backstep_nonl1_data.contribs_df.query('expected_editing_eff<=20')['sampletype'].tolist())
fwd_samples_set=set(fwd_samples)
rev_samples_set=set(rev_samples)
high_conc_fwd_samples=list(fwd_samples_set & high_conc_samples_set)
high_conc_rev_samples=list(rev_samples_set & high_conc_samples_set)
analyze_one_condition_data_custom(backstep_nonl1_data, high_conc_fwd_samples,  high_conc_rev_samples, low_samples, res_dir + '/analysis_result_highconc_backstep_nonl1_PtWAVE', 80, 100, 80, 100)
analyze_one_condition_data_custom(backstep_l1_data, high_conc_fwd_samples,  high_conc_rev_samples, low_samples, res_dir + '/analysis_result_highconc_backstep_l1_PtWAVE', 80, 100, 80, 100)
analyze_one_condition_data_custom(all_nonl1_data, high_conc_fwd_samples,  high_conc_rev_samples, low_samples, res_dir + '/analysis_result_highconc_all_nonl1_PtWAVE', 80, 100, 80, 100)
analyze_one_condition_data_custom(all_l1_data, high_conc_fwd_samples,  high_conc_rev_samples, low_samples, res_dir + '/analysis_result_highconc_all_l1_PtWAVE', 80, 100, 80, 100)
#!!!論文ポイント
# 高濃度下での傾向は実証できなかった、むしろnonl1のほうが優れている。またわずかではあるがallよりもbackstepのほうが精度が高い
low_conc_fwd_samples=list(fwd_samples_set & low_conc_samples_set)
low_conc_rev_samples=list(rev_samples_set & low_conc_samples_set)
analyze_one_condition_data_custom(backstep_nonl1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_backstep_nonl1_PtWAVE', 0, 20, 0, 20)
analyze_one_condition_data_custom(backstep_l1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_backstep_l1_PtWAVE', 0, 20, 0, 20)
analyze_one_condition_data_custom(all_nonl1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_all_nonl1_PtWAVE', 0, 20, 0, 20)
analyze_one_condition_data_custom(all_l1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_all_l1_PtWAVE', 0, 20, 0, 20)
#!!!論文ポイント
# 低濃度下での傾向はall/nonl1が優れている点は示すことができた。backstep/nonl1は相関としてはかなり微妙だが少なくともbackstep/l1よりはよい。
# l1は微量検出において悪影響を与えやすくbackstepと重なると相乗効果になるため、基本的には同時に選択すべきではない可能性がある
# この傾向は高濃度下でもみられるので、基本的に選択すべきではない組み合わせとなっている
### 通常DNA量 vs 低DNA量間比較
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), backstep_nonl1_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_backstep_nonl1')
comp_two_factor(random_nonl1_data.contribs_df.query('sampletype in ["5"]'), random_nonl1_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_random_nonl1')
comp_two_factor(all_nonl1_data.contribs_df.query('sampletype in ["5"]'), all_nonl1_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_all_nonl1')
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), backstep_nonl1_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_backstep_nonl1')
comp_two_factor(random_nonl1_data.contribs_df.query('sampletype in ["5"]'), random_nonl1_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_random_nonl1')
comp_two_factor(all_nonl1_data.contribs_df.query('sampletype in ["5"]'), all_nonl1_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_all_nonl1')
#!!!論文ポイント
# DNA量が少ないサンプル（低シグナルサンプル）では変異を過少に見積もる傾向がある。

### ツール間比較

# TIDEデフォルト解析データ
tide_default_data = SynDNAData(file_path='data_TIDEbatch_default.xlsx', external_analysis_mode='TIDE')
analyze_one_condition_data(tide_default_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_tide_default")
# TIDEロングデリーション設定解析データ
tide_longdel_setting_data = SynDNAData(file_path='data_TIDEbatch_size_50 decomposition_1_700 align_lb1.xlsx', external_analysis_mode='TIDE')
analyze_one_condition_data(tide_longdel_setting_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_tide_longdel_setting")
#!!!論文ポイント
# TIDEはもちろん-85欠失を検出できないが、効率推定の時点でほとんど相関がなく役に立たない。ラージデリーションが存在する場所で著しく効率を見誤る可能性が高い

# ICE解析データ
ice_data = SynDNAData(file_path='ICE/ICE_Results_x8mu9xn239uvg3qk/summary_x8mu9xn239uvg3qk.csv', external_analysis_mode='ICE')
analyze_one_condition_data(ice_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_ice")
#!!!論文ポイント
# TIDEと同様ラージデリーション検出では効率含めて大幅に値を見誤る

# DECODRバルク用設定解析データ
decodr_bulk_data = SynDNAData(file_path='data_DECODR_v3_0_bulk_standard_default.xlsx', external_analysis_mode='DECODR')
analyze_one_condition_data(decodr_bulk_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_decodr_bulk")
#!!!論文ポイント
# DECODRはバルク型だと100%野生型配列を変異型配列と判定することがある。これは致命的な点なのでPtWAVEでは適切に処理できたことも含めて論文で言及してもよい。
### 通常DNA量 vs 低DNA量間比較
comp_two_factor(decodr_bulk_data.contribs_df.query('sampletype in ["5"]'), decodr_bulk_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_decodr_bulk')
comp_two_factor(decodr_bulk_data.contribs_df.query('sampletype in ["5"]'), decodr_bulk_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_decodr_bulk')
### 50%濃度データのPtWAVE vs DECODR比較
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), decodr_bulk_data.contribs_df.query('sampletype in ["5"]'), 'editing_eff', 50, res_dir + '/innercomp_50per_ptwave_vs_decodr_bulk')

# DECODRクローン用設定解析データ
decodr_clonal_data = SynDNAData(file_path='data_DECODR_v3_0_clonal_standard_default.xlsx', external_analysis_mode='DECODR')
analyze_one_condition_data(decodr_clonal_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_decodr_clonal")
#!!!論文ポイント
# DECODRはクローンのほうがましな結果がでる。効率推定はクローン用解析がPtWAVEと同程度だが、-85欠失の検出精度ははるかに低い。これはDECODRのクローンモードでは-85欠失100%サンプルを-84欠失が100%と誤判定した点が大きい。
#### 高濃度変異サンプルの精度
analyze_one_condition_data_custom(decodr_clonal_data, high_conc_fwd_samples,  high_conc_rev_samples, low_samples, res_dir + '/analysis_result_highconc_decodr_clonal', 80, 100, 80, 100)
#!!!論文ポイント
# 高濃度変異サンプルでは0%検出があるため明らかに評価値が悪くなっている。グラフにするまでもないので決定係数・相関係数情報だけ記載する形でよいと思う。
#### 低濃度変異サンプルの精度
analyze_one_condition_data_custom(decodr_clonal_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_decodr_clonal', 0, 20, 0, 20)
#!!!論文ポイント
# 低濃度変異サンプルの検出のおいてはbackstep_nonl1には勝るが、all_nonl1よりも精度は悪い。この3条件でのデータもちゃんと示しておく。
### 通常DNA量 vs 低DNA量間比較
comp_two_factor(decodr_clonal_data.contribs_df.query('sampletype in ["5"]'), decodr_clonal_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_decodr_clonal')
comp_two_factor(decodr_clonal_data.contribs_df.query('sampletype in ["5"]'), decodr_clonal_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_decodr_clonal')
#!!!論文ポイント
# バルクの場合含めて、forwardの解析はできなかった。そのため低濃度データの解析の安定性ではPtWAVEに劣る。
### 50%濃度データのPtWAVE vs DECODR比較
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), decodr_clonal_data.contribs_df.query('sampletype in ["5"]'), 'editing_eff', 50, res_dir + '/innercomp_50per_ptwave_vs_decodr_clonal')
#!!!論文ポイント
# バルクの場合含めて、PtWAVEよりも低くインデル値を見積もる傾向がある
#!!!論文ポイント
# 総括としてDECODRは全体的に精度がやや劣っており、特に編集配列が高濃度下ではその傾向がより明瞭であった。
# また中濃度においては効率を低く見積もる傾向があり、低DNA量化での安定性も比較的低い。このことは本来あるはずのインデルを見逃すことに繋がっている可能性がある。
