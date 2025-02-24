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

# all/L1正則化なし/150bp欠失/50bp挿入まで
seventyfiveINSfifty_all_nonl1_patt = 'PtWAVE_*_on_seventyfiveINSfifty_all_[0-9]*-[1-3].ab1'
# 正規表現確認
for name in glob.glob(seventyfiveINSfifty_all_nonl1_patt):
  print(name)
# 問題なし
seventyfiveINSfifty_all_nonl1_data = SynDNAData(patt=seventyfiveINSfifty_all_nonl1_patt, max_deletion=150)
analyze_one_condition_data(seventyfiveINSfifty_all_nonl1_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_seventyfiveINSfifty_all_nonl1_PtWAVE")

# backstep_nonl1_dataで、BIC値が二極化しているので分けて解析も行ってみる。
backstep_nonl1_data.contribs_df.query('BIC<10000')['-85']
backstep_nonl1_data.contribs_df.query('BIC<10000')['other_indels']
backstep_nonl1_data.contribs_df.query('BIC>10000')['-85']
backstep_nonl1_data.contribs_df.query('BIC>10000')['other_indels']

### PtWAVEパラメータ間解析
### 手法間比較
comp_accuracies(all_nonl1_data.contribs_df, random_nonl1_data.contribs_df, backstep_nonl1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_nonl1_100del3ins_comp')
comp_accuracies(all_l1_data.contribs_df, random_l1_data.contribs_df, backstep_l1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_l1_100del3ins_comp')
comp_accuracies(seventyfive_all_nonl1_data.contribs_df, seventyfive_random_nonl1_data.contribs_df, seventyfive_backstep_nonl1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_seventyfive_nonl1_100del3ins_comp')
comp_accuracies(ten_all_nonl1_data.contribs_df, ten_random_nonl1_data.contribs_df, ten_backstep_nonl1_data.contribs_df, fwd_samples + rev_samples, res_dir + '/innercomp_ten_nonl1_100del3ins_comp')
### L1 vs nonL1間比較
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
low_conc_fwd_samples=list(fwd_samples_set & low_conc_samples_set)
low_conc_rev_samples=list(rev_samples_set & low_conc_samples_set)
analyze_one_condition_data_custom(backstep_nonl1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_backstep_nonl1_PtWAVE', 0, 20, 0, 20)
analyze_one_condition_data_custom(backstep_l1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_backstep_l1_PtWAVE', 0, 20, 0, 20)
analyze_one_condition_data_custom(all_nonl1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_all_nonl1_PtWAVE', 0, 20, 0, 20)
analyze_one_condition_data_custom(all_l1_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_all_l1_PtWAVE', 0, 20, 0, 20)
### 通常DNA量 vs 低DNA量間比較
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), backstep_nonl1_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_backstep_nonl1')
comp_two_factor(random_nonl1_data.contribs_df.query('sampletype in ["5"]'), random_nonl1_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_random_nonl1')
comp_two_factor(all_nonl1_data.contribs_df.query('sampletype in ["5"]'), all_nonl1_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_all_nonl1')
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), backstep_nonl1_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_backstep_nonl1')
comp_two_factor(random_nonl1_data.contribs_df.query('sampletype in ["5"]'), random_nonl1_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_random_nonl1')
comp_two_factor(all_nonl1_data.contribs_df.query('sampletype in ["5"]'), all_nonl1_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_all_nonl1')

### ツール間比較

# TIDEデフォルト解析データ
tide_default_data = SynDNAData(file_path='data_TIDEbatch_default.xlsx', external_analysis_mode='TIDE')
analyze_one_condition_data(tide_default_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_tide_default")
# TIDEロングデリーション設定解析データ
tide_longdel_setting_data = SynDNAData(file_path='data_TIDEbatch_size_50 decomposition_1_700 align_lb1.xlsx', external_analysis_mode='TIDE')
analyze_one_condition_data(tide_longdel_setting_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_tide_longdel_setting")

# ICE解析データ
ice_data = SynDNAData(file_path='ICE/ICE_Results_x8mu9xn239uvg3qk/summary_x8mu9xn239uvg3qk.csv', external_analysis_mode='ICE')
analyze_one_condition_data(ice_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_ice")

# DECODRバルク用設定解析データ
decodr_bulk_data = SynDNAData(file_path='data_DECODR_v3_0_bulk_standard_default.xlsx', external_analysis_mode='DECODR')
analyze_one_condition_data(decodr_bulk_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_decodr_bulk")

### 通常DNA量 vs 低DNA量間比較
comp_two_factor(decodr_bulk_data.contribs_df.query('sampletype in ["5"]'), decodr_bulk_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_decodr_bulk')
comp_two_factor(decodr_bulk_data.contribs_df.query('sampletype in ["5"]'), decodr_bulk_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_decodr_bulk')
### 50%濃度データのPtWAVE vs DECODR比較
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), decodr_bulk_data.contribs_df.query('sampletype in ["5"]'), 'editing_eff', 50, res_dir + '/innercomp_50per_ptwave_vs_decodr_bulk')

# DECODRクローン用設定解析データ
decodr_clonal_data = SynDNAData(file_path='data_DECODR_v3_0_clonal_standard_default.xlsx', external_analysis_mode='DECODR')
analyze_one_condition_data(decodr_clonal_data, fwd_samples, rev_samples, low_samples, res_dir + "/analysis_result_decodr_clonal")

#### 高濃度変異サンプルの精度
analyze_one_condition_data_custom(decodr_clonal_data, high_conc_fwd_samples,  high_conc_rev_samples, low_samples, res_dir + '/analysis_result_highconc_decodr_clonal', 80, 100, 80, 100)

#### 低濃度変異サンプルの精度
analyze_one_condition_data_custom(decodr_clonal_data, low_conc_fwd_samples,  low_conc_rev_samples, low_samples, res_dir + '/analysis_result_lowconc_decodr_clonal', 0, 20, 0, 20)

### 通常DNA量 vs 低DNA量間比較
comp_two_factor(decodr_clonal_data.contribs_df.query('sampletype in ["5"]'), decodr_clonal_data.contribs_df.query('sampletype in ["22"]'), 'editing_eff', 50, res_dir + '/innercomp_fwd_conccomp_decodr_clonal')
comp_two_factor(decodr_clonal_data.contribs_df.query('sampletype in ["5"]'), decodr_clonal_data.contribs_df.query('sampletype in ["23"]'), 'editing_eff', 50, res_dir + '/innercomp_rev_conccomp_decodr_clonal')

### 50%濃度データのPtWAVE vs DECODR比較
comp_two_factor(backstep_nonl1_data.contribs_df.query('sampletype in ["5"]'), decodr_clonal_data.contribs_df.query('sampletype in ["5"]'), 'editing_eff', 50, res_dir + '/innercomp_50per_ptwave_vs_decodr_clonal')
