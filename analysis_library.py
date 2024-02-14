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
from scipy.stats import pearsonr
from scipy.stats import ttest_1samp
import seaborn as sns
import datetime

class SynDNAData:
  contribs_df=None
  indel_ratios_stats=None
  max_deletion=None
  sample_exp_ref = {
    '1':0,
    '2':5,
    '3':10,
    '4':20,
    '5':50,
    '6':80,
    '7':90,
    '8':95,
    '9':100,
    '10':0,
    '11':0,
    '12':5,
    '13':80,
    '14':90,
    '15':95,
    '16':100,
    '17':0,
    '18':99,
    '19':99,
    '20':0,
    '21':0,
    '22':50,
    '23':50
  } 
  def __init__(self, patt, max_deletion):
    self.max_deletion = max_deletion
    recordable_range = [str(i) for i in list(range(-self.max_deletion, 51, 1))]
    self.contribs_df = pd.DataFrame(columns=recordable_range)
    for name in glob.glob(patt):
      name_nonext = os.path.splitext(name)[0]# ファイル名解析
      underscores = [i for i, char in enumerate(name_nonext) if char == '_']
      samplename = name_nonext[underscores[-1]:][1:]
      json_fn = os.path.join(os.getcwd(), name, "indel.json")# ファイルパス確認
      dash = [i for i, char in enumerate(samplename) if char == '-']# サンプル番号とレプリケート番号を抽出
      sampletype = samplename[:dash[-1]]
      samplerep = samplename[dash[-1]:][1:]
      exist_json_fn = os.path.isfile(json_fn)
      if exist_json_fn:
        with open(json_fn) as f:
          json_di = json.load(f)# jsonをデータフレームにコピー
          record_df = pd.DataFrame([json_di["editing_outcomes"]], index=[samplename])
          record_df['editing_eff'] = json_di["editing_eff"]
          record_df['r_sq'] = json_di["r_sq"]
          record_df['BIC'] = json_di["BIC"]
          record_df['sampletype'] = str(sampletype)
          record_df['samplerep'] = str(samplerep)
          record_df['expected_editing_eff'] = self.sample_exp_ref[str(sampletype)]
          print(record_df)
          self.contribs_df = pd.concat([self.contribs_df, record_df], ignore_index=False)# contribs_dfに情報をまとめる
          self.contribs_df = self.contribs_df.fillna(0.0)
      else:
        # record_df = pd.DataFrame(columns=recordable_range,index=[samplename])
        # record_df['sampletype'] = str(sampletype)
        # record_df['samplerep'] = str(samplerep)
        # record_df['expected_editing_eff'] = self.sample_exp_ref[str(sampletype)]
        # record_df['is_complete'] = False
        # print(record_df)
        continue
    self._Classify_Indels()
  def __init__(self, file_path, external_analysis_mode):
    if external_analysis_mode=='TIDE' or external_analysis_mode=='DECODR':
      excel_df = pd.read_excel(file_path)
      contribs_cols = [str(i) for i in excel_df['shiftrange'].to_list()]
      excel_rows = excel_df.columns.to_list()
      excel_rows.pop(0)
      contribs_rows = [i.replace('percentage ', '') for i in excel_rows]
      self.contribs_df = excel_df.set_axis(contribs_cols, axis=0).T.drop('shiftrange', axis=0).set_axis(contribs_rows, axis=0).drop('batch_analysis', axis=1)
      self.contribs_df['sampletype'] = [i.split('-')[0] for i in contribs_rows]
      self.contribs_df['samplerep'] = [i.split('-')[1] for i in contribs_rows]
      self.contribs_df['expected_editing_eff'] = [self.sample_exp_ref[i] for i in self.contribs_df['sampletype'].to_list()]
      self.max_deletion = 100
      print(self.contribs_df)
      self._Classify_Indels()
    elif external_analysis_mode=='ICE':
      csv_df = pd.read_csv(file_path)
      fltr_csv_df = csv_df[csv_df['ICE'].notna()]
      self.max_deletion = 100
      recordable_range = [str(i) for i in list(range(-self.max_deletion, 51, 1))]
      self.contribs_df = pd.DataFrame(columns=recordable_range)
      for id, row in fltr_csv_df.iterrows():
        indels_df = pd.DataFrame(json.loads(row['Indels'].replace('\'', '\"')), index=[0])
        record_df = pd.DataFrame(columns=recordable_range)
        record_df = pd.concat([record_df, indels_df], ignore_index=False)# contribs_dfに情報をまとめる
        record_df = record_df.fillna(0.0)
        sample_name = row['Label'].replace('s', '')
        record_df = record_df.set_axis([sample_name], axis=0)
        record_df['sampletype'] = sample_name.split('-')[0]
        record_df['samplerep'] = sample_name.split('-')[1]
        record_df['expected_editing_eff'] = [sample_exp_ref[i] for i in record_df['sampletype'].to_list()]
        record_df['editing_eff'] = row['ICE']
        record_df['r_sq'] = row['R Squared']
        record_df['BIC'] = None
        self.contribs_df = pd.concat([self.contribs_df, record_df], ignore_index=False)
      print(self.contribs_df)
      self._Classify_Indels()
  def _Classify_Indels(self):
    print("Classify indels...")
    other_indels_range = [str(i) for i in list(range(-self.max_deletion, -85, 1))] + [str(i) for i in list(range(-84, 0, 1))] + [str(i) for i in list(range(1, 51, 1))]
    print(other_indels_range)
    self.contribs_df['other_indels'] = self.contribs_df[other_indels_range].sum(axis=1)
    print(self.contribs_df)
    print("Calculate statistics...")
    print("Mean...")
    contribs_mean_df = self.contribs_df.groupby('sampletype').mean()
    print("SEM...")
    contribs_sem_df = self.contribs_df.groupby('sampletype').sem()
    print("Make Indel calss table...")
    indel_class_mean = contribs_mean_df[['-85','0', 'other_indels', 'expected_editing_eff']]
    indel_class_mean.rename(columns={'-85': '-85_mean', '0': '0_mean', 'other_indels': 'other_indels_mean'}, inplace=True)
    indel_class_sem = contribs_sem_df[['-85','0', 'other_indels']]
    indel_class_sem.rename(columns={'-85': '-85_sem', '0': '0_sem', 'other_indels': 'other_indels_sem'}, inplace=True)
    indel_class_mean['expected_editing_eff'] = round(contribs_mean_df['expected_editing_eff'])
    self.indel_class_stats = pd.concat([indel_class_mean, indel_class_sem], ignore_index=False, axis=1)
    print("Calculate statistics for analysis summary...")
    print("Make analysis summary table...")
    analysis_summary_mean = contribs_mean_df[['editing_eff','r_sq', 'BIC', 'expected_editing_eff']]
    analysis_summary_mean.rename(columns={'editing_eff': 'editing_eff_mean', 'r_sq': 'r_sq_mean', 'BIC': 'BIC_indels_mean'}, inplace=True)
    analysis_summary_sem = contribs_sem_df[['editing_eff','r_sq', 'BIC']]
    analysis_summary_sem.rename(columns={'editing_eff': 'editing_eff_sem', 'r_sq': 'r_sq_sem', 'BIC': 'BIC_indels_sem'}, inplace=True)
    analysis_summary_mean['expected_editing_eff'] = round(analysis_summary_mean['expected_editing_eff'])
    self.analysis_summary = pd.concat([analysis_summary_mean, analysis_summary_sem], ignore_index=False, axis=1)
  def SelectData(self, selection_list, indel_class_stats_path, analysis_summary_path):
    selected_indel_class_stats = self.indel_class_stats.query('sampletype in '+ str(selection_list))
    selected_indel_class_stats = selected_indel_class_stats.reindex(index=selection_list)
    selected_analysis_summary = self.analysis_summary.query('sampletype in '+ str(selection_list))
    selected_analysis_summary = selected_analysis_summary.reindex(index=selection_list)
    selected_indel_class_stats.to_excel(indel_class_stats_path, sheet_name='result')
    selected_analysis_summary.to_excel(analysis_summary_path, sheet_name='result')
    return (selected_indel_class_stats, selected_analysis_summary)

def analyze_one_condition_data(SynDNAData_instance, fwd_samples, rev_samples, low_samples, savedir):
  os.makedirs(savedir)
  fwd_indel_class_stats_df, fwd_analysis_summary_df = SynDNAData_instance.SelectData(fwd_samples, os.path.join(savedir, 'fwd_indel_class_stats.xlsx'), os.path.join(savedir, 'fwd_analysis_summary.xlsx'))
  rev_indel_class_stats_df, rev_analysis_summary_df = SynDNAData_instance.SelectData(rev_samples, os.path.join(savedir, 'rev_indel_class_stats.xlsx'), os.path.join(savedir, 'rev_analysis_summary.xlsx'))
  rev_indel_class_stats_df, rev_analysis_summary_df = SynDNAData_instance.SelectData(low_samples, os.path.join(savedir, 'low_indel_class_stats.xlsx'), os.path.join(savedir, 'low_analysis_summary.xlsx'))

  def actual_expected_comp_plot(specific_indel, editing_eff, expected_editing_eff, saveprefix):
    def specific_indel_plot(x, y, saveprefix):
      plt.figure(figsize=(5,5))
      plt.xlim(-5, 105)
      plt.xticks(np.arange(0,110,10))
      plt.ylim(-5, 105)
      plt.yticks(np.arange(0,110,10))
      specific_indel_reg = LinearRegression()
      specific_indel_avalue = x
      specific_indel_evalue = y
      specific_indel_reg.fit(specific_indel_avalue, specific_indel_evalue)
      plt.plot(specific_indel_avalue, specific_indel_evalue, 'rs', markersize=4)
      plt.plot(specific_indel_avalue, specific_indel_reg.predict(specific_indel_avalue), 'k-')
      r2_value = r2_score(specific_indel_avalue, specific_indel_evalue)
      correlation = pearsonr(specific_indel_avalue['-85'].to_list(), specific_indel_evalue['expected_editing_eff'].to_list())
      r_value = correlation[0]
      p_value = correlation[1]
      plt.savefig(saveprefix + '_specific_indel_r2=' + str(r2_value) + 'R=' + str(r_value) + 'p=' + str(p_value) + '.tiff', dpi=350)
      plt.clf()
    
    def editing_eff_plot(x, y, saveprefix):
      plt.figure(figsize=(5,5))
      plt.xlim(-5, 105)
      plt.xticks(np.arange(0,110,10))
      plt.ylim(-5, 105)
      plt.yticks(np.arange(0,110,10))
      editing_eff_reg = LinearRegression()
      editing_eff_avalue = x
      editing_eff_evalue = y
      editing_eff_reg.fit(editing_eff_avalue, editing_eff_evalue)
      plt.plot(editing_eff_avalue, editing_eff_evalue, 'bs', markersize=4)
      plt.plot(editing_eff_avalue, editing_eff_reg.predict(editing_eff_avalue), 'k-')
      r2_value = r2_score(editing_eff_avalue, editing_eff_evalue)
      correlation = pearsonr(editing_eff_avalue['editing_eff'].to_list(), editing_eff_evalue['expected_editing_eff'].to_list())
      r_value = correlation[0]
      p_value = correlation[1]
      plt.savefig(saveprefix + '_editing_eff_r2=' + str(r2_value) + 'R=' + str(r_value) + 'p=' + str(p_value) + '.tiff', dpi=350)
      plt.clf()
    specific_indel_plot(specific_indel, expected_editing_eff, saveprefix)
    editing_eff_plot(editing_eff, expected_editing_eff, saveprefix)

  def sample_indel_plot(contribs_df, savedir, max_deletion):
    indel_df = contribs_df[[str(i) for i in list(range(-max_deletion, 51, 1))]]
    for index, row in indel_df.iterrows():
      plt.figure(figsize=(23,5))
      plt.ylim(0, 100)
      row.plot(kind='bar', color='#D96704')
      plt.title(index)
      plt.savefig(savedir + index + '_indel_barplot', dpi=350)
      plt.clf()
  
  # 相関グラフをつくる
  normal_pmol_result_df = SynDNAData_instance.contribs_df.query('sampletype in '+ str(fwd_samples + rev_samples))
  actual_expected_comp_plot(normal_pmol_result_df[['-85']], normal_pmol_result_df[['editing_eff']], normal_pmol_result_df[['expected_editing_eff']], os.path.join(savedir, 'plot'))
  # インデル分布図を作成する
  sample_indel_plot(SynDNAData_instance.contribs_df, savedir + '/', SynDNAData_instance.max_deletion)

def analyze_one_condition_data_custom(SynDNAData_instance, fwd_samples, rev_samples, low_samples, savedir, xmin, xmax, ymin, ymax):
  os.makedirs(savedir)
  fwd_indel_class_stats_df, fwd_analysis_summary_df = SynDNAData_instance.SelectData(fwd_samples, os.path.join(savedir, 'fwd_indel_class_stats.xlsx'), os.path.join(savedir, 'fwd_analysis_summary.xlsx'))
  rev_indel_class_stats_df, rev_analysis_summary_df = SynDNAData_instance.SelectData(rev_samples, os.path.join(savedir, 'rev_indel_class_stats.xlsx'), os.path.join(savedir, 'rev_analysis_summary.xlsx'))
  rev_indel_class_stats_df, rev_analysis_summary_df = SynDNAData_instance.SelectData(low_samples, os.path.join(savedir, 'low_indel_class_stats.xlsx'), os.path.join(savedir, 'low_analysis_summary.xlsx'))
  def actual_expected_comp_plot(specific_indel, editing_eff, expected_editing_eff, saveprefix):
    def specific_indel_plot(x, y, saveprefix):
      plt.figure(figsize=(5,5))
      plt.xlim(xmin - 5, xmax + 5)
      plt.xticks(np.arange(xmin, xmax + 10,10))
      plt.ylim(ymin - 5, ymax + 5)
      plt.yticks(np.arange(ymin, ymax + 10,10))
      specific_indel_reg = LinearRegression()
      specific_indel_avalue = x
      specific_indel_evalue = y
      specific_indel_reg.fit(specific_indel_avalue, specific_indel_evalue)
      plt.plot(specific_indel_avalue, specific_indel_evalue, 'rs', markersize=4)
      plt.plot(specific_indel_avalue, specific_indel_reg.predict(specific_indel_avalue), 'k-')
      r2_value = r2_score(specific_indel_avalue, specific_indel_evalue)
      plt.savefig(saveprefix + '_specific_indel_r2=' + str(r2_value) + '.tiff', dpi=350)
      plt.clf()
    def editing_eff_plot(x, y, saveprefix):
      plt.figure(figsize=(5,5))
      plt.xlim(xmin - 5, xmax + 5)
      plt.xticks(np.arange(xmin, xmax + 10,10))
      plt.ylim(ymin - 5, ymax + 5)
      plt.yticks(np.arange(ymin, ymax + 10,10))
      editing_eff_reg = LinearRegression()
      editing_eff_avalue = x
      editing_eff_evalue = y
      editing_eff_reg.fit(editing_eff_avalue, editing_eff_evalue)
      plt.plot(editing_eff_avalue, editing_eff_evalue, 'bs', markersize=4)
      plt.plot(editing_eff_avalue, editing_eff_reg.predict(editing_eff_avalue), 'k-')
      r2_value = r2_score(editing_eff_avalue, editing_eff_evalue)
      plt.savefig(saveprefix + '_editing_eff_r2=' + str(r2_value) + '.tiff', dpi=350)
      plt.clf()
    specific_indel_plot(specific_indel, expected_editing_eff, saveprefix)
    editing_eff_plot(editing_eff, expected_editing_eff, saveprefix)
  def sample_indel_plot(contribs_df, savedir, max_deletion):
    indel_df = contribs_df[[str(i) for i in list(range(-max_deletion, 51, 1))]]
    for index, row in indel_df.iterrows():
      plt.figure(figsize=(23,5))
      plt.ylim(0, 100)
      row.plot(kind='bar', color='#D96704')
      plt.title(index)
      plt.savefig(savedir + index + '_indel_barplot', dpi=350)
      plt.clf()
  normal_pmol_result_df = SynDNAData_instance.contribs_df.query('sampletype in '+ str(fwd_samples + rev_samples))# 相関グラフをつくる
  actual_expected_comp_plot(normal_pmol_result_df[['-85']], normal_pmol_result_df[['editing_eff']], normal_pmol_result_df[['expected_editing_eff']], os.path.join(savedir, 'plot'))
  # sample_indel_plot(SynDNAData_instance.contribs_df, savedir + '/', SynDNAData_instance.max_deletion)# インデル分布図を作成する <-今の関数の用途では不要

def comp_accuracies(all_contribs_df, random_contribs_df, backstep_contribs_df, sample_list, savedir):
  os.makedirs(savedir)
  all_df = all_contribs_df.query('sampletype in '+ str(sample_list))
  random_df = random_contribs_df.query('sampletype in '+ str(sample_list))
  backstep_df = backstep_contribs_df.query('sampletype in '+ str(sample_list))
  def plot_three_box(all_arr, random_arr, backstep_arr, y_label, ymax, savename):
    plt.figure(figsize=(5,10))
    sns.set()
    sns.set_style('whitegrid')
    sns.set_palette('Set3')
    df = pd.DataFrame({
      'All': all_arr.to_list(),
      'Random': random_arr.to_list(),
      'Backstep': backstep_arr.to_list()
    })
    df_melt = pd.melt(df, var_name='method', value_name=y_label)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sns.boxplot(x='method', y=y_label, data=df_melt, showfliers=False, ax=ax)
    sns.stripplot(x='method', y=y_label, data=df_melt, jitter=True, color='black', ax=ax, linewidth=1, alpha=.3)
    plt.ylim(0, ymax)
    plt.savefig(os.path.join(savedir, savename + '.tiff'), dpi=350)
    plt.clf()
  total_arr = all_df['BIC'].to_list() + random_df['BIC'].to_list() + backstep_df['BIC'].to_list()
  plot_three_box(all_df['r_sq'], random_df['r_sq'], backstep_df['r_sq'], 'R2', 1.0, 'Accurary')
  plot_three_box(all_df['BIC'], random_df['BIC'], backstep_df['BIC'], 'BIC', max(total_arr)+1000, 'Unstability')
  stat_df = pd.DataFrame({
    'all_r2_mean': all_df[['r_sq']].mean(),
    'all_r2_sem': all_df[['r_sq']].sem(),
    'all_bic_mean': all_df[['BIC']].mean(),
    'all_bic_sem': all_df[['BIC']].sem(),
    'random_r2_mean': random_df[['r_sq']].mean(),
    'random_r2_sem': random_df[['r_sq']].sem(),
    'random_bic_mean': random_df[['BIC']].mean(),
    'random_bic_sem': random_df[['BIC']].sem(),
    'backstep_r2_mean':backstep_df[['r_sq']].mean(),
    'backstep_r2_sem': backstep_df[['r_sq']].sem(),
    'backstep_bic_mean': backstep_df[['BIC']].mean(),
    'backstep_bic_sem': backstep_df[['BIC']].sem(),
    'r2_all_random_pvalue': stats.ttest_ind(all_df[['r_sq']], random_df[['r_sq']], equal_var=False).pvalue[0],
    'r2_random_backstep_pvalue': stats.ttest_ind(random_df[['r_sq']], backstep_df[['r_sq']], equal_var=False).pvalue[0],
    'r2_backstep_all_pvalue': stats.ttest_ind(backstep_df[['r_sq']], all_df[['r_sq']], equal_var=False).pvalue[0],
    'bic_all_random_pvalue': stats.ttest_ind(all_df[['BIC']], random_df[['BIC']], equal_var=False).pvalue[0],
    'bic_random_backstep_pvalue': stats.ttest_ind(random_df[['BIC']], backstep_df[['BIC']], equal_var=False).pvalue[0],
    'bic_backstep_all_pvalue': stats.ttest_ind(backstep_df[['BIC']], all_df[['BIC']], equal_var=False).pvalue[0],
    'r2_all_random_nonpara_pvalue': stats.wilcoxon(all_df['r_sq'], random_df['r_sq'], alternative="two-sided").pvalue,
    'r2_random_backstep_nonpara_pvalue': stats.wilcoxon(random_df['r_sq'], backstep_df['r_sq'], alternative="two-sided").pvalue,
    'r2_backstep_all_nonpara_pvalue': stats.wilcoxon(backstep_df['r_sq'], all_df['r_sq'], alternative="two-sided").pvalue,
    'bic_all_random_nonpara_pvalue': stats.wilcoxon(all_df['BIC'], random_df['BIC'], alternative="two-sided").pvalue,
    'bic_random_backstep_nonpara_pvalue': stats.wilcoxon(random_df['BIC'], backstep_df['BIC'], alternative="two-sided").pvalue,
    'bic_backstep_all_nonpara_pvalue': stats.wilcoxon(backstep_df['BIC'], all_df['BIC'], alternative="two-sided").pvalue
  })
  stat_df.to_excel(os.path.join(savedir, 'stat.xlsx'), sheet_name='result')

def comp_two_factor(A_contribs_df, B_contribs_df, target_label, expected_value, savedir):
  os.makedirs(savedir)
  def plot_two_box(A_arr, B_arr, y_label, ymax, expected_value, savename):
    plt.figure(figsize=(5,10))
    sns.set()
    sns.set_style('whitegrid')
    sns.set_palette('Set3')
    df = pd.DataFrame({
      'A': A_arr.to_list(),
      'B': B_arr.to_list(),
    })
    df_melt = pd.melt(df, var_name='method', value_name=y_label)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    sns.boxplot(x='method', y=y_label, data=df_melt, showfliers=False, ax=ax)
    sns.stripplot(x='method', y=y_label, data=df_melt, jitter=True, color='black', ax=ax, linewidth=1, alpha=.3)
    plt.axhline(expected_value, color='red')
    plt.ylim(0, ymax)
    plt.savefig(os.path.join(savedir, savename + '.tiff'), dpi=350)
    plt.clf()
  total_arr = A_contribs_df[target_label].to_list() + B_contribs_df[target_label].to_list()
  plot_two_box(A_contribs_df[target_label], B_contribs_df[target_label], target_label, max(total_arr)*1.1, expected_value, target_label)
  stat_df = pd.DataFrame({
    'A_mean': [A_contribs_df[target_label].mean()],
    'A_sem': [A_contribs_df[target_label].sem()],
    'B_mean': [B_contribs_df[target_label].mean()],
    'B_sem': [B_contribs_df[target_label].sem()],
    'A_one_ttest_pvalue':[ttest_1samp(A_contribs_df[target_label], popmean=expected_value).pvalue],
    'B_one_ttest_pvalue':[ttest_1samp(B_contribs_df[target_label], popmean=expected_value).pvalue],
    'welch_ttest_pvalue': [stats.ttest_ind(A_contribs_df[target_label].to_list(), B_contribs_df[target_label].to_list(), equal_var=False).pvalue],
    'nonpara_pvalue': [stats.wilcoxon(A_contribs_df[target_label].to_list(), B_contribs_df[target_label].to_list(), alternative="two-sided").pvalue]
  })
  stat_df.to_excel(os.path.join(savedir, 'stat.xlsx'), sheet_name='result')