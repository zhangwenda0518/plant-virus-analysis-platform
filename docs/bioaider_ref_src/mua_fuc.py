# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: mutation_tools/mua_fuc.py
"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2024/7/18 19:57

"""
import matplotlib.pyplot as plt
import numpy as np
from custom_libraries.public_functions import make_dir, aa_kind_dis, identify_gap_unknow, mutation_index

def piecewise_read(input_file_path, max_num_seq, run_nums):
    """
    分段读取数据到内存进行计算
    :param input_file_path: 输入序列名称
    :param max_num_seq: 每次最大读取的序列条数
    :param run_nums: 第几次读取
    :return:m +1, seq_tab, run_flage，特别注意seq_tab里还带了序列下标，在其他分析可能用不到
    """
    seq_tab = []
    seq_dic = {}
    seq_name_list = []
    seq_name = ""
    with open(input_file_path, "r", encoding="utf-8") as gen_file:
        for line in gen_file:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line.strip(">")
                seq_name_list.append(seq_name)
                seq_dic[seq_name] = []
            elif line != "":
                seq_dic[seq_name].append(line)

        n_start = max_num_seq * (run_nums - 1)
        n_end = max_num_seq * run_nums
        m = len(seq_name_list) - 1
        for n in range(len(seq_name_list)):
            if n_start <= n < n_end:
                each_seqname = seq_name_list[n]
                seq = "".join(seq_dic[each_seqname]).upper()
                seq_tab.append([n + 1, each_seqname, seq])

        if m < n_end:
            run_flage = False
        else:
            run_flage = True
    return (
     m + 1, seq_tab, run_flage)


def identify_degenerate_bases(base):
    """
    判断是否含有简并碱基
    :param base: 核苷酸序列片段
    :return: False or True
    """
    degenerate_bases = {
     'W', 'Y', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V'}
    degenerate = False
    if base in degenerate_bases:
        degenerate = True
    return degenerate


def qujianpuanduan(shizhi_str, group_list):
    """
    突变频率分组_区间判断
    :param shizhi_str: 突变频率数值
    :param group_list: 突变频率分组
    :return: 该突变频率数值位于的分组
    """
    qujian = "Others"
    shizhi = int(shizhi_str)
    for group in group_list:
        group = group.strip()
        k1, k2 = group.split("\t")
        k1 = int(k1)
        k2 = int(k2)
        if k1 < shizhi <= k2:
            qujian = "(" + str(k1) + "," + str(k2) + "]"
            break

    return qujian


def link_mutation_run1(seq_name, seq_lenth, seq_contain, qury_seq):
    """
     适用于核苷酸和氨基酸关联突变分析
    :param seq_name: 序列名称
    :param seq_lenth: 序列长度
    :param seq_contain: 序列碱基
    :param qury_seq: 参考序列碱基
    :return:
    """
    seq_mutation_list = []
    for site_number in range(seq_lenth):
        analy_site = seq_contain[site_number[:site_number + 1]]
        qury_seq_site = qury_seq[site_number[:site_number + 1]]
        if analy_site != qury_seq_site:
            seq_mutation_list.append(qury_seq_site + str(site_number + 1) + analy_site)

    mutation_str = " + ".join(seq_mutation_list)
    if mutation_str != "":
        return [
         seq_name, mutation_str]


def link_mutation_run2(seq_name, seq_contain, qury_seq, cycles):
    """
    适用于密码子关联突变分析
    :param seq_name:
    :param seq_contain:
    :param qury_seq:
    :param cycles:
    :return:
    """
    seq_mutation_list = []
    for codon_number in range(cycles):
        analy_site = seq_contain[(3 * codon_number)[:3 * (codon_number + 1)]]
        qury_seq_site = qury_seq[(3 * codon_number)[:3 * (codon_number + 1)]]
        if analy_site != qury_seq_site:
            seq_mutation_list.append(qury_seq_site + str(codon_number + 1) + analy_site)

    mutation_str = " + ".join(seq_mutation_list)
    if mutation_str != "":
        return [
         seq_name, mutation_str]


def single_mutation_nt(seq_name, qury_seq, seq_contain, seq_lenth):
    """
    核苷酸突变分析（位点扫描模式下）
    :param seq_name: 该条序列名称
    :param qury_seq: 参考序列碱基
    :param seq_contain: 该条序列的碱基
    :param seq_lenth: 该条序列的长度
    :return:
    """
    site_mution = []
    site_mution_JX = []
    for site in range(seq_lenth):
        qury_seq_site = qury_seq[site[:site + 1]]
        s = seq_contain[site[:site + 1]]
        if s != qury_seq_site:
            site_mution.append([qury_seq_site, str(site + 1), s])
            if s != "-" and qury_seq_site != "-" and identify_degenerate_bases(s) == False and identify_degenerate_bases(qury_seq_site) == False:
                site_mution_JX.append([qury_seq_site, str(site + 1), s])

    return [
     seq_name, site_mution, site_mution_JX]


def single_mutation_codon(qury_seq, seq_name, seq, codon_count, codon_table, translate):
    """
    密码子突变分析（位点扫描模式下）
    :param qury_seq: 参考序列密码子
    :param seq_name: 该条序列名称
    :param seq: 该条序列密码子
    :param codon_count: 密码子总数
    :param codon_table: 翻译表
    :param translate: 翻译方法
    :return:
    """
    mutation_codon_list = []
    for site in range(codon_count):
        qury_seq_conding = qury_seq[(3 * site)[:3 * (site + 1)]]
        s = seq[(3 * site)[:3 * (site + 1)]]
        if s != qury_seq_conding:
            site = site + 1
            qury_seq_aa = translate(qury_seq_conding, codon_table)
            s_aa = translate(s, codon_table)
            if qury_seq_conding == "---":
                mutation = "\t".join([qury_seq_conding, str(site), s, "Insertion", "-", s_aa])
            else:
                if s == "---":
                    mutation = "\t".join([qury_seq_conding, str(site), s, "Deletion", qury_seq_aa, "-"])
                else:
                    if qury_seq_aa == s_aa:
                        mutation = "\t".join([qury_seq_conding, str(site), s, "Synonymous", qury_seq_aa, s_aa])
                    else:
                        if s_aa == "?":
                            mutation = "\t".join([qury_seq_conding, str(site), s, "Unknown", qury_seq_aa, s_aa])
                        else:
                            if s_aa == "*":
                                mutation = "\t".join([qury_seq_conding, str(site), s, "Termination", qury_seq_aa, s_aa])
                            else:
                                mutation = "\t".join([qury_seq_conding, str(site), s, "Nonsynonymous", qury_seq_aa, s_aa])
            mutation_codon_list.append(mutation)

    return [seq_name, mutation_codon_list]


def save_codon_log(savefile_path, qury_seq_name, mutation_results_list):
    """
    密码子突变分析（位点扫描模式下）生成日志文件
    :param savefile_path: 日志文件输出路径
    :param qury_seq_name: 该条序列的名称
    :param mutation_results_list: 突变位点列表
    :return:
    """
    with open((savefile_path + "_run_logs.txt"), "a", encoding="utf-8") as log_txt:
        log_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "Sequence name" + "\t" + " Site (codon)" + "\t" + "Index of nt in codon" + "\t" + "mutation type" + "\t" + "mutation" + "\n")
        for each_seq in mutation_results_list:
            seq_result = each_seq.get()
            seq_name = seq_result[0]
            mutation_list = seq_result[1]
            if mutation_list != []:
                for each_mutation in mutation_list:
                    mutation_codon_list = each_mutation.split("\t")
                    mutation_nt_location = mutation_index(mutation_codon_list[0], mutation_codon_list[2])
                    mutation_str = "".join([mutation_codon_list[0], " to ", mutation_codon_list[2],
                     " (", mutation_codon_list[-2], mutation_codon_list[1], mutation_codon_list[-1], ")"])
                    logs = "\t".join([seq_name, mutation_codon_list[1],
                     mutation_nt_location,
                     mutation_codon_list[3],
                     mutation_str])
                    log_txt.write(logs + "\n")


def single_mutation_aa(qury_seq, seq_lenth, seq_name, seq_contain):
    """
    氨基酸突变分析（位点扫描模式下）
    :param qury_seq: 参考序列氨基酸
    :param seq_lenth: 该条序列长度（对齐后）
    :param seq_name: 该条序列名称
    :param seq_contain: 该条序列氨基酸
    :return:
    """
    aa_mutation = []
    aa_mutation_JX = []
    for site in range(seq_lenth):
        qury_seq_site = qury_seq[site[:site + 1]]
        s = seq_contain[site[:site + 1]]
        if s != qury_seq_site:
            aa_mutation.append([qury_seq_site, str(site + 1), s, aa_kind_dis(qury_seq_site, s)])
            if identify_gap_unknow(s) == False and identify_gap_unknow(qury_seq_site) == False:
                aa_mutation_JX.append([qury_seq_site, str(site + 1), s, aa_kind_dis(qury_seq_site, s)])

    return [
     seq_name, aa_mutation, aa_mutation_JX]


def nuc_aa_plot(site_fru, frequency_group_tab, savefile_path):
    """
    核苷酸或氨基酸突变分析绘制频率分布图
    :param site_fru: 突变频率表
    :param frequency_group_tab: 突变频率分组表
    :param savefile_path: 输出文件路径
    :return:
    """
    dd = {}
    for every_site in site_fru:
        furequence_type = qujianpuanduan(site_fru[every_site], frequency_group_tab)
        key = furequence_type
        if key not in dd:
            dd[key] = 1
        else:
            dd[key] += 1

    for each_group in frequency_group_tab:
        k1, k2 = each_group.split("\t")
        new_key = "(" + str(k1) + "," + str(k2) + "]"
        flage = 0
        for keys in dd:
            if keys == new_key:
                flage = 1
                break

        if flage == 0:
            dd[new_key] = 0

    fretab_save_path = savefile_path + "_frequency_groups_of_substitution.txt"
    woqu = open(fretab_save_path, "w", encoding="utf-8")
    woqu.write("Frequency_groups\tCount\n")
    fretab_file = []
    for furequence_interval in dd:
        woqu.write(furequence_interval + "\t" + str(dd[furequence_interval]) + "\n")
        fretab_file.append(furequence_interval + "\t" + str(dd[furequence_interval]))

    woqu.close()
    group_infoma = []
    for line in frequency_group_tab:
        jj1, jj2 = line.split("\t")
        linn = "(" + str(jj1) + "," + str(jj2) + "]"
        group_infoma.append(linn)

    fre_conut = []
    for x in group_infoma:
        for line in fretab_file:
            fruence, numbers = line.strip().split("\t")
            if fruence == x:
                fre_conut.append(int(numbers))
                break

    label_count = len(group_infoma)
    bar_width = 0.4
    plt.bar((np.arange(label_count)), fre_conut, align="center",
      color="#00C5CD",
      alpha=0.8,
      width=bar_width,
      edgecolor="k")
    plt.xlabel("Group of substitution frequency")
    plt.ylabel("Substitution sites count")
    plt.title("Frequency distribution of substitution sites")
    plt.xticks((np.arange(label_count)), group_infoma, rotation=45)
    for x, y in enumerate(fre_conut):
        plt.text(x, (y + 0.1), ("%s" % round(y, 1)), ha="center")

    plt.savefig((savefile_path + "_substitution_frequency_distribution.png"),
      format="png",
      dpi=1200)
    plt.savefig((savefile_path + "_substitution_frequency_distribution.pdf"),
      format="pdf")
    plt.close()


def codon_plot(saves_dir, syn_nonsy_site, frequency_group_tab, out_prefix):
    """
    基于密码子突变分析绘制频率分布图
    :param saves_dir: 输出目录
    :param syn_nonsy_site: 同义替换和非同义替换位点
    :param frequency_group_tab: 频率分组表
    :param out_prefix: 输出文件的前缀
    :return:
    """
    base_statistics_draw_dir = saves_dir + "/" + "draw"
    make_dir(base_statistics_draw_dir)
    syn_nonsy_draw = open((base_statistics_draw_dir + "/nt_sites_for_draw.csv"), "w",
      encoding="utf-8")
    syn_nonsy_draw.write("Codon,Base index in codon,Nucleotide site,Type,Substitution frequency\n")
    for key in syn_nonsy_site:
        uf1, uf2 = key.split("\t")
        codon_site, nt_index = uf1.split("-")
        nt_num_index = 0
        if nt_index == "1":
            nt_num_index = 3 * (int(codon_site) - 1) + 1
        else:
            if nt_index == "2":
                nt_num_index = 3 * (int(codon_site) - 1) + 2
            else:
                if nt_index == "3":
                    nt_num_index = 3 * int(codon_site)
        syn_nonsy_draw.write(codon_site + "," + nt_index + "," + str(nt_num_index) + "," + uf2 + "," + str(syn_nonsy_site[key]) + "\n")

    syn_nonsy_draw.close()
    shizhi = open((base_statistics_draw_dir + "/nt_sites_for_draw.csv"), "r",
      encoding="utf-8")
    fu = shizhi.read().strip().split("\n")
    shizhi.close()
    dd = {}
    fretab_save_path = base_statistics_draw_dir + "/" + "Frequency_group_of_sys_and_non-sys_substitution.txt"
    woqu = open(fretab_save_path, "w", encoding="utf-8")
    woqu.write("Frequency_group\tSubstitution_type\tCount\n")
    for n in range(1, len(fu)):
        u1, u2, u3, u4, u5 = fu[n].split(",")
        furequence_type = qujianpuanduan(u5, frequency_group_tab)
        key = furequence_type + "\t" + u4
        if key not in dd:
            dd[key] = 1
        else:
            dd[key] += 1

    alltype = {}
    group_infos = []
    for gro in frequency_group_tab:
        gro = gro.strip()
        kk1, kk2 = gro.split("\t")
        qujianj = "(" + str(kk1) + "," + str(kk2) + "]"
        group_infos.append(qujianj)
        key1 = qujianj + "\t" + "Synonymous"
        key2 = qujianj + "\t" + "Nonsynonymous"
        alltype[key1] = 0
        alltype[key2] = 0

    fretab_file = []
    for key in alltype:
        flage = 0
        for key2 in dd:
            if key2 == key:
                woqu.write(key2 + "\t" + str(dd[key2]) + "\n")
                fretab_file.append(key2 + "\t" + str(dd[key2]))
                flage = 1
                break

        if flage == 0:
            woqu.write(key + "\t" + str(0) + "\n")
            fretab_file.append(key + "\t" + str(0))

    woqu.close()
    nonsy_fre_conut = []
    sy_fre_conut = []
    for x in group_infos:
        for line in fretab_file:
            fruence, sus_type, numbers = line.strip().split("\t")
            if fruence == x:
                if sus_type == "Nonsynonymous":
                    nonsy_fre_conut.append(int(numbers))
                if sus_type == "Synonymous":
                    sy_fre_conut.append(int(numbers))

    label_count = len(group_infos)
    bar_width = 0.4
    plt.bar((np.arange(label_count)), nonsy_fre_conut, label="Nonsynonymous",
      color="#EE5C42",
      alpha=0.8,
      width=bar_width,
      edgecolor="k")
    plt.bar((np.arange(label_count) + bar_width), sy_fre_conut,
      label="Synonymous",
      color="#00C5CD",
      alpha=0.8,
      width=bar_width,
      edgecolor="k")
    plt.xlabel("Group of substitution frequency")
    plt.ylabel("Substitution nt sites count")
    plt.title("Frequency distribution of substitution sites")
    plt.xticks((np.arange(label_count) + bar_width / 2), group_infos,
      rotation=45)
    for x, y in enumerate(nonsy_fre_conut):
        plt.text(x, y + 0.1, "%s" % y)

    for x, y in enumerate(sy_fre_conut):
        plt.text(x + bar_width, y + 0.1, "%s" % y)

    plt.legend()
    plt.savefig((base_statistics_draw_dir + "/" + out_prefix + "_substitution_frequency_distribution.png"),
      format="png",
      dpi=1200)
    plt.savefig((base_statistics_draw_dir + "/" + out_prefix + "_substitution_frequency_distribution.pdf"),
      format="pdf")
    plt.close()


def lollipop_data(lollipop_infro):
    lollipop_single_nt = {}
    for each_line in lollipop_infro:
        codon_site, base_inx, kinds, mutation_infro, seqcount = each_line
        mutations = mutation_infro.split(" (")[0]
        if base_inx.find(",") != -1:
            base_inx_list = base_inx.split(",")
            for base in base_inx_list:
                nt_num_index = 0
                if base == "1":
                    nt_num_index = 3 * (int(codon_site) - 1) + 1
                else:
                    if base == "2":
                        nt_num_index = 3 * (int(codon_site) - 1) + 2
                    else:
                        if base == "3":
                            nt_num_index = 3 * int(codon_site)
                codon_change = mutations.split(" to ")
                mutation_record = "".join([
                 codon_change[0][int(base) - 1],
                 str(nt_num_index),
                 codon_change[1][int(base) - 1]])
                key = "\t".join([str(nt_num_index), kinds, mutation_record])
                if key in lollipop_single_nt:
                    lollipop_single_nt[key] += int(seqcount)
                else:
                    lollipop_single_nt[key] = int(seqcount)

        else:
            nt_num_index = 0
            if base_inx == "1":
                nt_num_index = 3 * (int(codon_site) - 1) + 1
            else:
                if base_inx == "2":
                    nt_num_index = 3 * (int(codon_site) - 1) + 2
                else:
                    if base_inx == "3":
                        nt_num_index = 3 * int(codon_site)
            codon_change = mutations.split(" to ")
            mutation_record = "".join([
             codon_change[0][int(base_inx) - 1],
             str(nt_num_index),
             codon_change[1][int(base_inx) - 1]])
            key = "\t".join([str(nt_num_index), kinds, mutation_record])
            if key in lollipop_single_nt:
                lollipop_single_nt[key] += int(seqcount)
            else:
                lollipop_single_nt[key] = int(seqcount)

    return lollipop_single_nt
