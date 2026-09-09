# -*- coding: utf-8 -*-
"""Web 层共用小工具（多 blueprint 共享，自 app.py 拆出）。"""
import os

from vp.utils import check_path
from vp.config import DIRS


def _safe_sample(sample):
    import re
    return re.sub(r'[^A-Za-z0-9_\-.]', '_', str(sample))
