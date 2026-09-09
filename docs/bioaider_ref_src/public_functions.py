# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: custom_libraries/public_functions.py
"""
Author: Zhou Zhi-Jian

"""
import sys, os, platform, re, config
app_dir = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))
if platform.system().lower() == "windows":
    app_dir = app_dir.replace("\\", "/")

def load_update_plugin(events):
    """
    加载插件,更新插件路径
    :param events:
    :return:
    """
    plugins_list_file_path = app_dir + "/Configurations/plug_set.init"
    if events == "load":
        with open(plugins_list_file_path, "r", encoding="utf-8") as plugins_list_file:
            for line in plugins_list_file:
                if not line.startswith("#"):
                    line_list = line.split("\t")
                    config.plugins_list.append([line_list[0], line_list[1], line_list[2], line_list[3], line_list[4].strip()])

    else:
        if events == "update":
            plugins_list_file = open(plugins_list_file_path, "w", encoding="utf-8")
            plugins_list_file.write("#Name\tStatus\tDescription\tOperation\tLocal path\n")
            for each_plugin in config.plugins_list:
                plugins_list_file.write("\t".join(each_plugin) + "\n")


def set_font_size(operation_type, new_font):
    """
    获取字体，字号
    :param operation_type:
    :param new_font:
    :return:
    """
    setfile_path = app_dir + "/Configurations/others.init"
    if operation_type.upper() == "OPEN":
        app_font = "Calibri,12"
        with open(setfile_path, "r", encoding="utf-8") as font_size_file:
            for line in font_size_file:
                if line.startswith("app_font"):
                    app_font = re.findall("{(.*?)}", line)[0]

        return app_font
    if operation_type.upper() == "WRITE":
        with open(setfile_path, "r", encoding="utf-8") as font_size_file:
            file_contain = font_size_file.read()
        with open(setfile_path, "w", encoding="utf-8") as font_size_file:
            new_contain = re.sub("app_font = {(.*?)}", "app_font = {" + new_font + "}", file_contain, 1)
            font_size_file.write(new_contain)


def make_dir(input_dir):
    """
    :param input_dir: 创建文件夹，如果该文件夹存在就忽视
    :return:
    """
    if os.path.isdir(input_dir) == False:
        os.makedirs(input_dir)


def resolve_file_path(file_path):
    """
    解析文件路径的目录，后缀
    :param file_path: 文件路径
    :return:
    """
    input_data_dir = os.path.dirname(file_path)
    inputfile_name = file_path.split("/")[-1]
    file_farmat = inputfile_name.split(".")[-1]
    out_prefix = inputfile_name.replace(file_farmat, "").replace(".", "_")
    return (
     input_data_dir, out_prefix)


def get_all_path(open_dir_path):
    """
    获取某一文件夹里所有文件的路径 (包括子文件夹)
    :param open_dir_path: 文件夹
    :return: 所有文件的路径列表
    """
    rootdir = open_dir_path
    path_list = []
    lists = os.listdir(rootdir)
    for i in range(0, len(lists)):
        com_path = os.path.join(rootdir, lists[i])
        com_path = com_path.replace("\\", "/")
        if os.path.isfile(com_path):
            path_list.append(com_path)
        if os.path.isdir(com_path):
            path_list.extend(get_all_path(com_path))

    return path_list


def change_app_theme(app_stylesheet_init_path, slect_theme):
    """
    更改软件主题
    :param app_stylesheet_init_path:
    :param slect_theme:
    :return:
    """
    app_stylesheet_init = open(app_stylesheet_init_path, "w",
      encoding="utf-8")
    app_stylesheet_init.write("stylesheet:" + slect_theme)
    app_stylesheet_init.close()
    return


def isnum(n):
    """
    判断是否为浮点数
    :param n: 数字或者字符
    :return: True or False
    """
    try:
        t = float(n)
        return True
    except:
        return False


def fasta_all_to_one(filepath):
    """
    多行fasta序列转一行标准格式(制作标准的fasta格式序列）
    :param filepath: fasta序列文件路径
    :return: 标准的fasta格式序列字符串
    """
    file_contains = []
    with open(filepath, "r", encoding="utf-8") as gen_file:
        for line in gen_file:
            line = line.strip()
            if line.startswith(">"):
                file_contains.append("\n" + line + "\n")
            else:
                file_contains.append(line)

    return "".join(file_contains).strip()


def fast_read_fasta(input_file_path):
    """
    读取序列为列表
    :param input_file_path:
    :return:
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

        for each_seqname in seq_name_list:
            seq = "".join(seq_dic[each_seqname]).upper()
            seq_tab.append([each_seqname, seq])

    return seq_tab


def read_fasta_dic(input_file_path):
    """
    超快读取fasta序列，储存为字典,1.8Gb宿主数据用时18秒。8.8万条新冠N蛋白序列用时0.5秒。
    :param input_file_path: fasta序列文件路径
    :return: 返回一个字典（key为序列名称，value为序列）
    """
    seq_dic = {}
    seq_name = ""
    with open(input_file_path, "r", encoding="utf-8") as gen_file:
        for line in gen_file:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line
                seq_dic[seq_name] = []
            elif line != "":
                seq_dic[seq_name].append(line.upper())

        for key, value in seq_dic.items():
            seq_dic[key] = "".join(value)

    return seq_dic


def fas_to_fas(seq):
    """
    fasta转fasta格式
    :param seq: fasta格式序列
    :return: 标准fasta格式序列
    """
    a = seq.split("\n")
    contains = []
    for line in a:
        line = line.strip()
        if line != "":
            if line.startswith(">"):
                contains.append("\n" + line + "\n")
            else:
                contains.append(line)

    return "".join(contains).strip()


def fas_to_nex(seq, datatype):
    """
    fasta转nexus格式(连续)
    :param seq: 为标准的fasta格式的字符串
    :param datatype: DNA 或者 PROTEIN
    :return: 连续nexus格式的序列
    """
    f2 = seq.split("\n")
    seqname_list = []
    seq = []
    for line in f2:
        line = line.strip()
        if line != "":
            if line.startswith(">"):
                seqname_list.append(line.replace(">", "").replace(" ", "_"))
            else:
                seq.append(line)

    ntax = len(seqname_list)
    maxlenname = max((len(i) for i in seqname_list))
    nchar = max((len(i) for i in seq))
    nexus_title = "#NEXUS\n\nBegin data;\n"
    para1 = "\tDimensions ntax=" + str(ntax) + " nchar=" + str(nchar) + ";" + "\n"
    para2 = "\tFormat datatype=" + datatype + " gap=- missing=?;" + "\n"
    para3 = "\tMatrix\n"
    end_flage = "\t;\nEnd;"
    body = []
    for n in range(len(seqname_list)):
        seqtaxon = seqname_list[n] + " " * (maxlenname - len(seqname_list[n]))
        body.append(seqtaxon + "    " + seq[n] + "\n")

    nexus_seq = nexus_title + para1 + para2 + para3 + "".join(body) + end_flage
    return nexus_seq


def fas_to_phy(seq):
    """
    fasta转phylip
    :param seq: 标准的fasta格式的字符串
    :return: phylip的序列
    """
    f2 = seq.split("\n")
    Seqname = []
    seq = []
    for line in f2:
        line = line.strip()
        if line != "":
            if line.startswith(">"):
                Seqname.append(line.replace(">", "").replace(" ", "_"))
            else:
                seq.append(line)

    ntax = len(Seqname)
    maxlenname = max((len(i) for i in Seqname))
    nchar = max((len(i) for i in seq))
    phylip_title = str(ntax) + " " + str(nchar) + "\n"
    body = []
    for n in range(len(Seqname)):
        seqtaxon = Seqname[n] + " " * (maxlenname - len(Seqname[n]))
        body.append(seqtaxon + "    " + seq[n] + "\n")

    phylip_seq = phylip_title + "".join(body)
    return phylip_seq


def fas_to_pml(seq):
    """
    fasta转paml格式
    :param seq: 标准的fasta格式的字符串
    :return: paml格式序列
    """
    f2 = seq.split("\n")
    Seqname = []
    seq = []
    for line in f2:
        line = line.strip()
        if line != "":
            if line.startswith(">"):
                Seqname.append(line.replace(">", "").replace(" ", "_"))
            else:
                seq.append(line)

    ntax = len(Seqname)
    maxlenname = max((len(i) for i in Seqname))
    nchar = max((len(i) for i in seq))
    paml_title = str(ntax) + "\t" + str(nchar) + "\n"
    body = []
    for n in range(len(Seqname)):
        seqtaxon = Seqname[n] + " " * (maxlenname - len(Seqname[n]))
        body.append(seqtaxon + "\n" + seq[n] + "\n")

    paml_seq = paml_title + "".join(body)
    return paml_seq


def nex_to_fas(seq):
    """
    nexus转fasta
    :param seq: 连续nexus格式的序列
    :return: 标准的fasta格式序列
    """
    seq = re.findall("Matrix(.*?)End;", seq.replace("\n", "&"), re.I)[0]
    seq = seq.replace("&", "\n")
    seq = seq.replace(";", "")
    seq = seq.strip()
    f2 = seq.split("\n")
    fasta_seq = []
    for line in f2:
        line = line.replace("\t", "")
        line = line.strip()
        fasta_seq.append(">" + line.split(" ")[0] + "\n" + line.split(" ")[-1] + "\n")

    return "".join(fasta_seq)


def phy_to_fas(seq):
    """
    phylip转fasta
    :param seq: phylip格式的序列字符串
    :return: 标准fasta格式序列
    """
    cons = seq.split("\n")
    fasta_seq = []
    for n in range(1, len(cons)):
        cons[n] = cons[n].strip()
        if cons[n] != "":
            fasta_seq.append("".join([">", cons[n].split(" ")[0], "\n", cons[n].split(" ")[-1], "\n"]))

    return "".join(fasta_seq)


def pml_to_fas(seq):
    """
    paml转fasta
    :param seq:  paml格式的序列字符串
    :return:  标准的fasta格式序列
    """
    cons = seq.split("\n")
    cons.remove(cons[0])
    k = len(cons)
    fasta_seq = []
    for i in range(int(k / 2)):
        fasta_seq.append("".join([">", cons[2 * i], "\n", cons[2 * i + 1].upper()]))

    return "\n".join(fasta_seq)


def aa_kind_dis(a1, a2):
    """
    判断非同义突变对应的氨基酸极性和带电性变化
    :param a1: 参考氨基酸
    :param a2: 突变后的氨基酸
    :return: 性质改变
    """
    AA_inftab = {
     'G': '"(polar,none-charge)"', 'S': '"(polar,none-charge)"', 
     'T': '"(polar,none-charge)"', 'C': '"(polar,none-charge)"', 
     'Q': '"(polar,none-charge)"', 'N': '"(polar,none-charge)"', 
     'Y': '"(polar,none-charge)"', 'K': '"(polar,positive-charge)"', 
     'R': '"(polar,positive-charge)"', 'H': '"(polar,positive-charge)"', 
     'D': '"(polar,negative-charge)"', 'E': '"(polar,negative-charge)"', 
     'A': '"(non-polar,none-charge)"', 'V': '"(non-polar,none-charge)"', 
     'L': '"(non-polar,none-charge)"', 'I': '"(non-polar,none-charge)"', 
     'P': '"(non-polar,none-charge)"', 'F': '"(non-polar,none-charge)"', 
     'W': '"(non-polar,none-charge)"', 'M': '"(non-polar,none-charge)"', 
     '*': '"(None,None)"', '-': '"(None,None)"', 'X': '"(None,None)"'}
    if a1 == "?" or a2 == "?":
        aa_ty = "Unknown\tUnknown"
    else:
        b1 = AA_inftab[a1]
        b2 = AA_inftab[a2]
        if b1.split(",")[0] == b2.split(",")[0]:
            if b1.split(",")[1] == b2.split(",")[1]:
                aa_ty = "No\tNo"
            else:
                aa_ty = "".join(['Yes', '\t', b1, ' to ', b2])
        else:
            aa_ty = "".join(['Yes', '\t', b1, ' to ', b2])
    return aa_ty


def mutation_index(str1, str2):
    """
    定位密码子中的突变碱基位点
    :param str1: 参考密码子
    :param str2: 突变密码子
    :return: 碱基位点
    """
    nt_index = []
    for n in range(3):
        if str1[n[:n + 1]] != str2[n[:n + 1]]:
            nt_index.append(str(n + 1))

    return ",".join(nt_index)


def identify_gap_unknow(aa):
    """
    判断蛋白质序列中是否含有gap、*、X和？
    :param aa: 氨基酸序列片段
    :return: False or True
    """
    gap_unknow = {
     "-", "*", "?", "X"}
    identify_result = False
    if aa in gap_unknow:
        identify_result = True
    return identify_result


def delete_gap(seq):
    """
    删除fasta序列中的gap
    :param seq: 普通fasta格式的序列
    :return: 不含gap的标准fasta格式序列
    """
    a = seq.split("\n")
    contains = []
    for line in a:
        line = line.strip()
        if line.startswith(">"):
            contains.append("\n" + line + "\n")
        else:
            contains.append(line.replace("-", ""))

    return "".join(contains).strip()


def cut_str(seqstr, k1, k2):
    """
    定义获取子字符串函数
    :param seqstr: 字符串
    :param k1: 最小子串长度
    :param k2: 最大子串长度
    :return: 子串列表
    """
    results = []
    for x in range(k1 - 1, k2):
        for i in range(len(seqstr) - x):
            results.append(seqstr[i[:i + x + 1]])

    cuts = list(set(results))
    cuts.sort(key=(results.index))
    return cuts


def repeat_search(qure_str, min_len, max_len):
    strs = cut_str(qure_str, min_len, max_len)
    result_list = []
    for x in strs:
        start = [x.start() + 1 for x in re.finditer(x, qure_str)]
        counts = len(start)
        propo = counts * len(x) / len(qure_str)
        sites_list = []
        for k in start:
            sites_list.append(str(k))

        sites = ",".join(sites_list)
        if sites.find(",") != -1:
            result_list.append([x, sites, str(propo)])

    return result_list


def single_find(re_pattern, context):
    """
    返回正则表达式匹配结果
    :param re_pattern: 匹配模式
    :param context: 用于匹配的文本
    :return: 匹配结果
    """
    results = "Unknown"
    find_list = re_pattern.findall(context)
    if find_list == []:
        results = results
    else:
        results = find_list[0].strip()
    return results


def date_to_float(date_str):
    """
    日期转浮点数
    :param date_str: 输入的日期，比如2020.10.28
    :param rules: 日期格式
    :return: 输出结果，浮点数
    """
    date_days1 = [
     31, 28, 31, 30, 
     31, 30, 31, 31, 
     30, 
     31, 30, 31]
    date_days2 = [
     31, 29, 31, 30, 
     31, 30, 31, 31, 
     30, 
     31, 30, 31]
    try:
        try:
            date_list = re.split("[./_-]", date_str)
            if len(date_list) == 3:
                month = int(date_list[1])
                if int(date_list[0]) % 4 == 0:
                    days = sum([date_days2[i] for i in range(0, month - 1)]) + int(date_list[2])
                    results = str(float(date_list[0]) + days / 366.00000001)
                else:
                    days = sum([date_days1[i] for i in range(0, month - 1)]) + int(date_list[2])
                    results = str(float(date_list[0]) + days / 365.00000001)
            else:
                if len(date_list) == 2:
                    date_list = re.split("[./_-]", date_str)
                    month = int(date_list[1])
                    if int(date_list[0]) % 4 == 0:
                        days = sum([date_days2[i] for i in range(0, month - 1)]) + 1
                        results = str(float(date_list[0]) + days / 366.00000001)
                    else:
                        days = sum([date_days1[i] for i in range(0, month - 1)]) + 1
                        results = str(float(date_list[0]) + days / 365.00000001)
                else:
                    results = date_str
        except:
            results = date_str

    finally:
        pass

    return results
