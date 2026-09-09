# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: custom_libraries/calculationlib.py
"""
Author: Zhou Zhi-Jian
Time: 2023/2/27 15:17

"""
import textdistance

def seqIndetity(str1, str2):
    """
    计算序列同一性
    :param str1: 字符串1
    :param str2: 字符串2
    :return: 同一性数值
    """
    K = int(len(str1) - len(str2))
    if K > 0:
        str2 = str2 + "-" * K
    else:
        if K < 0:
            str1 = str1 + "-" * -K
    seq_a1 = []
    seq_a2 = []
    for x1 in str1:
        seq_a1.append(x1)

    for x2 in str2:
        seq_a2.append(x2)

    pairgap_count = 0
    matching_conut = 0
    mismatch_conut = 0
    for n in range(len(seq_a1)):
        if seq_a1[n] == seq_a2[n]:
            if seq_a1[n] == "-":
                pairgap_count = pairgap_count + 1
            matching_conut = matching_conut + 1
        else:
            mismatch_conut = mismatch_conut + 1

    identity_str1_str2 = (matching_conut - pairgap_count) / (matching_conut - pairgap_count + mismatch_conut)
    return identity_str1_str2


def seqIndetity_Condense_gap(str1, str2):
    """
    计算序列同一性(考虑压缩gap的情况下)
    :param str1: 字符串1
    :param str2: 字符串2
    :return: 同一性数值
    """
    K = int(len(str1) - len(str2))
    if K > 0:
        str2 = str2 + "-" * K
    else:
        if K < 0:
            str1 = str1 + "-" * -K
    seq_a1 = []
    seq_a2 = []
    for x1 in str1:
        seq_a1.append(x1)

    for x2 in str2:
        seq_a2.append(x2)

    pairgap_count = 0
    matching_conut = 0
    mismatch_conut = 0
    inser = 0
    for n in range(len(seq_a1)):
        if seq_a1[n] == seq_a2[n]:
            if seq_a1[n] == "-":
                pairgap_count = pairgap_count + 1
            matching_conut = matching_conut + 1

    identity_str1_str2 = (matching_conut - pairgap_count) / (matching_conut - pairgap_count + (mismatch_conut - 2 * inser / 3))
    return identity_str1_str2


def similary_counter(way_tye, strn, strm):
    """
    计算序列相似性
    :param way_tye: 选择计算similary的算法
    :param strn: 字符串1（序列1）
    :param strm: 字符串2（序列2）
    :return: 相似性值
    """
    pair_similary = float(0)
    if way_tye == "Sequence Identity":
        pair_similary = seqIndetity(strn, strm)
    elif way_tye == "Hamming":
        pair_similary = textdistance.hamming.normalized_similarity(strn, strm)
    else:
        if way_tye == "Smith-Waterman":
            pair_similary = textdistance.smith_waterman.normalized_similarity(strn, strm)
        else:
            if way_tye == "Levenshtein":
                pair_similary = textdistance.levenshtein.normalized_similarity(strn, strm)
            else:
                if way_tye == "Damerau-Levenshtein":
                    pair_similary = textdistance.damerau_levenshtein.normalized_similarity(strn, strm)
                else:
                    if way_tye == "Needleman-Wunsch":
                        pair_similary = textdistance.needleman_wunsch.normalized_similarity(strn, strm)
    return pair_similary
