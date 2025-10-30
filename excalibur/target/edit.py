'''target task inputs and parameter overwrite'''

# Heritage code shame:
# pylint: disable=too-many-lines

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import logging

log = logging.getLogger(__name__)


# ------------------- ------------------------------------------------
# -- CREATE -- -------------------------------------------------------
def createversion():
    '''
    1.2.0: Added Spitzer targets
    1.3.0: Added JWST-NIRISS-NIS-CLEAR-G700XD
    1.4.0: Added more targets (242 total)
    1.5.0: New limb darkening coefficients for spitzer targest
    1.6.0: Added a few FEH* values and one Hmag; removed some redundant settings
    1.7.0: Added confirmed planets in Ariel target list
    1.7.1: WFC3 targets 2023
    1.8.0: JWST filters
    1.8.1: I.M HST run
    '''
    return dawgie.VERSION(1, 8, 1)


# (NOTES FROM ADDITIONAL OF 500-some NEW TARGETS. Oct 2024)

# EPIC 205950854 is K2-168
# removed it from here and from run scripts
# ah wait this is a really weird one.  b is K2-168 but c is EPIC
#  (the Archive and POMA both use this odd mixmatch)
# If we just use K2-168, do we get 'c'?  yes, we're good!
#  (download table names do not match the main table names)
#  (that also explains how TOI-784 is in main table, but HD 307842 in download table)
#
# LHS 475 is GJ 4102
# LHS 475 : LHS475
# removed it from here and from run scripts
#
# HD 3167 : HD3167   is K2-96
# removed it from here and from run scripts
#
# Kepler-460 is KIC 5437945
#  this is another weird one (like K2-168 above) with two star names
#  Kepler-460 c but KIC 5437945 b  yikes

# WD 1856+534 is WD 1856   that's a bad name. let's switch it
#  removed these older ones:
# WD 1856 : WD 1856+534
# WD 1856 : WD1856

# TIC 172900988 is TIC 172900988 Aa
#  that's not a great name.  let's just use that as an Exoplanet Archive alias
#
# GJ 12 is Gliese 12
#  SIMBAD doesn't even list that name. just use it as an Exoplanet Archive alias
#
# TOI-784 is HD 307842 in the Archive
#  keep the TOI name  (arg! the aliases are only updated by target.create)

# oof there's a batch that work for target() but crash on system()
# they are all listed as false positives on the Exoplanet Archive
#  (let Jenn Burt know about these)
# all removed from scripts
# Kepler-488
# Kepler-494
# Kepler-628
# Kepler-706
# Kepler-807
# 4/7/25: also Kepler-470 and K2-399 (based on Sept.2024 paper)


# ------------ -------------------------------------------------------
# -- TARGET LIST -- --------------------------------------------------
# FIRST COL HAS TO BE SOLVABLE BY
# -- Obsolete https://archive.stsci.edu/hst/
# https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
# SECOND COL [OPTIONAL] HAS TO BE 'Default Alias' RETURNED FROM
# https://exoplanetarchive.ipac.caltech.edu/index.html
# OR ALIAS KNOWN BY https://archive.stsci.edu/hst/
# MINOR CHANGE
def targetlist():
    '''
    55 Cnc :
    CoRoT-1 :
    CoRoT-2 :
    CoRoT-7 :
    GJ 1132 :
    GJ 1214 :
    GJ 1252 :
    GJ 3053 : LHS 1140
    GJ 3470 :
    GJ 436 :
    GJ 9827 :
    HAT-P-1 :
    HAT-P-11 :
    HAT-P-12 :
    HAT-P-13 :
    HAT-P-15 :
    HAT-P-16 :
    HAT-P-17 :
    HAT-P-18 :
    HAT-P-19 :
    HAT-P-2 :
    HAT-P-20 :
    HAT-P-22 :
    HAT-P-23 :
    HAT-P-26 :
    HAT-P-3 :
    HAT-P-30 :
    HAT-P-32 :
    HAT-P-33 :
    HAT-P-34 :
    HAT-P-38 :
    HAT-P-4 :
    HAT-P-40 :
    HAT-P-41 :
    HAT-P-5 :
    HAT-P-6 :
    HAT-P-7 :
    HAT-P-8 :
    HATS-28 :
    HATS-3 :
    HATS-65 :
    HATS-7 :
    HD 106315 :
    HD 149026 :
    HD 17156 :
    HD 189733 :
    HD 209458 :
    HD 213885 :
    HD 219134 :
    HD 219666 :
    HD 23472 :
    HD 97658 :
    HR 858 :
    K2-124 :
    K2-132 :
    K2-136 :
    K2-138 :
    K2-141 :
    K2-167 :
    K2-174 :
    K2-18 :
    K2-19 :
    K2-21 :
    K2-212 :
    K2-22 :
    K2-24 :
    K2-25 :
    K2-26 :
    K2-28 :
    K2-289 :
    K2-3 :
    K2-31 :
    K2-32 :
    K2-33 :
    K2-36 :
    K2-52 :
    K2-53 :
    K2-55 :
    K2-58 :
    K2-79 :
    K2-87 :
    K2-9 :
    K2-90 :
    K2-93 : HIP 41378
    K2-95 :
    K2-96 : HD 3167
    K2-97 :
    KELT-1 :
    KELT-11 :
    KELT-14 :
    KELT-16 :
    KELT-20 :
    KELT-3 :
    KELT-7 :
    KELT-9 :
    Kepler-10 :
    Kepler-11 :
    Kepler-102 :
    Kepler-104 :
    Kepler-1083 :
    Kepler-12 :
    Kepler-125 :
    Kepler-126 :
    Kepler-127 :
    Kepler-13 :
    Kepler-1339 :
    Kepler-138 :
    Kepler-14 :
    Kepler-1485 :
    Kepler-1492 :
    Kepler-156 :
    Kepler-1568 :
    Kepler-158 :
    Kepler-16 :
    Kepler-1625 :
    Kepler-1651 :
    Kepler-167 :
    Kepler-17 :
    Kepler-18 :
    Kepler-19 :
    Kepler-20 :
    Kepler-205 :
    Kepler-218 :
    Kepler-236 :
    Kepler-249 :
    Kepler-25 :
    Kepler-26 :
    Kepler-293 :
    Kepler-297 :
    Kepler-309 :
    Kepler-32 :
    Kepler-37 :
    Kepler-395 :
    Kepler-45 :
    Kepler-454 :
    Kepler-48 :
    Kepler-482 :
    Kepler-49 :
    Kepler-5 :
    Kepler-504 :
    Kepler-505 :
    Kepler-570 :
    Kepler-582 :
    Kepler-598 :
    Kepler-6 :
    Kepler-603 :
    Kepler-61 :
    Kepler-62 :
    Kepler-68 :
    Kepler-7 :
    Kepler-705 :
    Kepler-737 :
    Kepler-769 :
    Kepler-786 :
    Kepler-9 :
    Kepler-93 :
    Kepler-94 :
    LHS 3844 :
    OGLE-TR-056 : OGLE-TR-56
    OGLE-TR-10 :
    Qatar-1 :
    Qatar-2 :
    TOI-270 :
    TOI-700 :
    TOI-849 :
    TrES-1 :
    TrES-2 :
    TrES-3 :
    TRAPPIST-1 :
    WASP-1 :
    WASP-10 :
    WASP-100 :
    WASP-101 :
    WASP-103 :
    WASP-104 :
    WASP-107 :
    WASP-11 :
    WASP-12 :
    WASP-121 :
    WASP-127 :
    WASP-13 :
    WASP-131 :
    WASP-14 :
    WASP-140 :
    WASP-15 :
    WASP-16 :
    WASP-17 :
    WASP-18 :
    WASP-19 :
    WASP-2 :
    WASP-21 :
    WASP-24 :
    WASP-26 :
    WASP-28 :
    WASP-29 :
    WASP-3 :
    WASP-31 :
    WASP-32 :
    WASP-33 :
    WASP-34 :
    WASP-35 :
    WASP-36 :
    WASP-37 :
    WASP-38 :
    WASP-39 :
    WASP-4 :
    WASP-43 :
    WASP-46 :
    WASP-48 :
    WASP-49 :
    WASP-5 :
    WASP-50 :
    WASP-52 :
    WASP-6 :
    WASP-62 :
    WASP-63 :
    WASP-64 :
    WASP-65 :
    WASP-67 :
    WASP-69 :
    WASP-7 :
    WASP-72 :
    WASP-74 :
    WASP-75 :
    WASP-76 :
    WASP-77 : WASP-77 A
    WASP-78 :
    WASP-79 :
    WASP-8 :
    WASP-80 :
    WASP-87 :
    WASP-94 : WASP-94 A
    WASP-95 :
    WASP-96 :
    WASP-97 :
    XO-1 :
    XO-2 : XO-2 N
    XO-3 :
    XO-4 :
    XO-5 :
    AU Mic :
    CoRoT-5 :
    CoRoT-11 :
    CoRoT-19 :
    DS Tuc A :
    EPIC 211945201 :
    EPIC 246851721 :
    G 9-40 :
    GJ 3473 :
    GJ 357 :
    GJ 367 :
    GJ 3929 :
    GJ 486 :
    GPX-1 :
    HAT-P-14 :
    HAT-P-21 :
    HAT-P-24 :
    HAT-P-25 :
    HAT-P-27 :
    HAT-P-28 :
    HAT-P-29 :
    HAT-P-31 :
    HAT-P-35 :
    HAT-P-36 :
    HAT-P-37 :
    HAT-P-39 :
    HAT-P-42 :
    HAT-P-43 :
    HAT-P-44 :
    HAT-P-45 :
    HAT-P-46 :
    HAT-P-49 :
    HAT-P-50 :
    HAT-P-51 :
    HAT-P-52 :
    HAT-P-53 :
    HAT-P-54 :
    HAT-P-55 :
    HAT-P-56 :
    HAT-P-57 :
    HAT-P-58 :
    HAT-P-59 :
    HAT-P-60 :
    HAT-P-61 :
    HAT-P-62 :
    HAT-P-64 :
    HAT-P-65 :
    HAT-P-66 :
    HAT-P-67 :
    HAT-P-68 :
    HAT-P-69 :
    HAT-P-70 :
    HAT-P-9 :
    HATS-1 :
    HATS-11 :
    HATS-13 :
    HATS-18 :
    HATS-2 :
    HATS-23 :
    HATS-24 :
    HATS-25 :
    HATS-26 :
    HATS-27 :
    HATS-29 :
    HATS-30 :
    HATS-31 :
    HATS-33 :
    HATS-34 :
    HATS-35 :
    HATS-37 A :
    HATS-38 :
    HATS-39 :
    HATS-4 :
    HATS-40 :
    HATS-41 :
    HATS-42 :
    HATS-43 :
    HATS-46 :
    HATS-47 :
    HATS-48 A :
    HATS-5 :
    HATS-50 :
    HATS-51 :
    HATS-52 :
    HATS-53 :
    HATS-56 :
    HATS-57 :
    HATS-58 A :
    HATS-6 :
    HATS-60 :
    HATS-62 :
    HATS-64 :
    HATS-67 :
    HATS-68 :
    HATS-70 :
    HATS-72 :
    HATS-9 :
    HD 108236 :
    HD 110082 :
    HD 110113 :
    HD 136352 :
    HD 1397 :
    HD 152843 :
    HD 15337 :
    HD 183579 :
    HD 191939 :
    HD 202772 A :
    HD 207897 :
    HD 221416 :
    HD 2685 :
    HD 332231 :
    HD 5278 :
    HD 63433 :
    HD 63935 :
    HD 73583 :
    HD 86226 :
    HD 89345 :
    HIP 65 A :
    HIP 67522 :
    K2-107 :
    K2-116 :
    K2-121 :
    K2-129 :
    K2-139 :
    K2-140 :
    K2-155 :
    K2-198 :
    K2-222 :
    K2-232 :
    K2-237 :
    K2-238 :
    K2-239 :
    K2-260 :
    K2-261 :
    K2-266 :
    K2-280 :
    K2-284 :
    K2-287 :
    K2-29 :
    K2-295 :
    K2-329 :
    K2-333 :
    K2-334 :
    K2-34 :
    K2-353 :
    K2-39 :
    K2-403 :
    K2-405 :
    K2-406 :
    KELT-10 :
    KELT-12 :
    KELT-15 :
    KELT-17 :
    KELT-18 :
    KELT-19 A :
    KELT-2 A :
    KELT-21 :
    KELT-23 A :
    KELT-24 :
    KELT-4 A :
    KELT-6 :
    KELT-8 :
    KOI-13 :
    KOI-94 :
    KPS-1 :
    Kepler-105 :
    Kepler-108 :
    Kepler-1314 :
    Kepler-1513 :
    Kepler-33 :
    Kepler-396 :
    Kepler-42 :
    Kepler-435 :
    Kepler-444 :
    Kepler-447 :
    Kepler-450 :
    Kepler-468 :
    Kepler-489 :
    Kepler-76 :
    Kepler-79 :
    L 98-59 :
    LHS 1478 :
    LHS 1678 :
    LP 714-47 :
    LP 791-18 :
    LTT 1445 A :
    LTT 3780 :
    LTT 9779 :
    MASCARA-1 :
    MASCARA-4 :
    NGTS-10 :
    NGTS-11 :
    NGTS-12 :
    NGTS-13 :
    NGTS-2 :
    NGTS-5 :
    NGTS-6 :
    NGTS-8 :
    Qatar-10 :
    Qatar-4 :
    Qatar-5 :
    Qatar-6 :
    Qatar-7 :
    Qatar-8 :
    Qatar-9 :
    TIC 257060897 :
    TOI-1064 :
    TOI-1075 :
    TOI-1130 :
    TOI-1201 :
    TOI-122 :
    TOI-1227 :
    TOI-1231 :
    TOI-125 :
    TOI-1259 A :
    TOI-1260 :
    TOI-1266 :
    TOI-1268 :
    TOI-1296 :
    TOI-1298 :
    TOI-1333 :
    TOI-1411 :
    TOI-1431 :
    TOI-1442 :
    TOI-1478 :
    TOI-150 :
    TOI-1518 :
    TOI-157 :
    TOI-1601 :
    TOI-163 :
    TOI-1670 :
    TOI-1685 :
    TOI-169 :
    TOI-1693 :
    TOI-172 :
    TOI-1728 :
    TOI-1759 :
    TOI-178 :
    TOI-1789 :
    TOI-1807 :
    TOI-1842 :
    TOI-1860 :
    TOI-1899 :
    TOI-201 :
    TOI-2076 :
    TOI-2109 :
    TOI-216 :
    TOI-2260 :
    TOI-2337 :
    TOI-237 :
    TOI-2411 :
    TOI-2427 :
    TOI-257 :
    TOI-2669 :
    TOI-269 :
    TOI-3362 :
    TOI-421 :
    TOI-431 :
    TOI-4329 :
    TOI-451 :
    TOI-481 :
    TOI-500 :
    TOI-530 :
    TOI-540 :
    TOI-544 :
    TOI-559 :
    TOI-561 :
    TOI-564 :
    TOI-620 :
    TOI-628 :
    TOI-640 :
    TOI-674 :
    TOI-677 :
    TOI-776 :
    TOI-813 :
    TOI-824 :
    TOI-833 :
    TOI-837 :
    TOI-892 :
    TOI-905 :
    TOI-954 :
    TrES-4 :
    TrES-5 :
    V1298 Tau :
    WASP-105 :
    WASP-106 :
    WASP-110 :
    WASP-113 :
    WASP-114 :
    WASP-117 :
    WASP-118 :
    WASP-119 :
    WASP-120 :
    WASP-123 :
    WASP-124 :
    WASP-126 :
    WASP-132 :
    WASP-133 :
    WASP-135 :
    WASP-136 :
    WASP-138 :
    WASP-139 :
    WASP-141 :
    WASP-142 :
    WASP-145 A :
    WASP-147 :
    WASP-148 :
    WASP-151 :
    WASP-153 :
    WASP-156 :
    WASP-157 :
    WASP-158 :
    WASP-159 :
    WASP-160 B :
    WASP-161 :
    WASP-163 :
    WASP-164 :
    WASP-165 :
    WASP-166 :
    WASP-167 :
    WASP-168 :
    WASP-169 :
    WASP-170 :
    WASP-172 :
    WASP-173 A :
    WASP-174 :
    WASP-175 :
    WASP-176 :
    WASP-177 :
    WASP-178 :
    WASP-180 A :
    WASP-181 :
    WASP-182 :
    WASP-183 :
    WASP-184 :
    WASP-185 :
    WASP-186 :
    WASP-187 :
    WASP-189 :
    WASP-190 :
    WASP-192 :
    WASP-20 :
    WASP-22 :
    WASP-23 :
    WASP-25 :
    WASP-41 :
    WASP-42 :
    WASP-44 :
    WASP-45 :
    WASP-47 :
    WASP-53 :
    WASP-54 :
    WASP-55 :
    WASP-56 :
    WASP-57 :
    WASP-58 :
    WASP-61 :
    WASP-66 :
    WASP-68 :
    WASP-70 A :
    WASP-71 :
    WASP-73 :
    WASP-81 :
    WASP-82 :
    WASP-83 :
    WASP-84 :
    WASP-85 A :
    WASP-88 :
    WASP-89 :
    WASP-90 :
    WASP-91 :
    WASP-92 :
    WASP-93 :
    WASP-98 :
    WASP-99 :
    Wolf 503 :
    XO-6 :
    XO-7 :
    pi Men : HD 39091
    HATS-22 :
    K2-30 :
    Kepler-1308 :
    Qatar-3 :
    WASP-129 :
    WASP-144 :
    Kepler-51 :
    WD 1856+534 :
    GJ 4102 : LHS 475
    HD 80606 :
    GJ 4332 : L 168-9
    LTT 5972 : TOI-836
    GJ 3090 :
    GJ 806 :
    HD 109833 :
    HD 207496 :
    HD 260655 :
    HD 93963 A :
    HD 95338 :
    HIP 29442 : HD 42813
    HIP 94235 :
    HIP 97166 :
    K2-105 :
    K2-370 :
    K2-415 :
    K2-60 :
    Kepler-1656 :
    Kepler-63 :
    Kepler-96 :
    NGTS-14 A :
    TOI-1136 :
    TOI-1272 :
    TOI-1278 :
    TOI-1288 :
    TOI-132 :
    TOI-139 :
    TOI-1422 :
    TOI-1468 :
    TOI-1470 :
    TOI-1634 :
    TOI-1694 :
    TOI-1695 :
    TOI-1710 :
    TOI-1801 :
    TOI-181 :
    TOI-1853 :
    TOI-1859 :
    TOI-199 :
    TOI-2000 :
    TOI-2010 :
    TOI-2018 :
    TOI-2134 :
    TOI-2136 :
    TOI-220 :
    TOI-2364 :
    TOI-2443 :
    TOI-2445 :
    TOI-2459 :
    TOI-2498 :
    TOI-251 :
    TOI-277 :
    TOI-3082 :
    TOI-332 :
    TOI-3629 :
    TOI-3785 :
    TOI-4010 :
    TOI-444 :
    TOI-4479 :
    TOI-4600 :
    TOI-4641 :
    TOI-470 :
    TOI-5126 :
    TOI-532 :
    TOI-5344 :
    TOI-5398 :
    TOI-5678 :
    TOI-5704 :
    TOI-5803 :
    TOI-672 :
    TOI-712 :
    TOI-908 :
    TOI-913 :
    TOI-942 :
    TOI-969 :
    Wolf 327 :
    CoRoT-3 :
    CoRoT-36 :
    Gaia-1 :
    Gaia-2 :
    HATS-10 :
    HATS-12 :
    HATS-45 :
    HATS-55 :
    HATS-61 :
    HD 118203 :
    HD 15906 :
    HD 21749 : GJ 143
    HD 235088 :
    HD 28109 :
    HIP 113103 :
    HIP 116454 :
    HIP 9618 :
    K2-233 :
    K2-240 :
    K2-285 :
    K2-321 :
    K2-344 :
    K2-348 :
    K2-417 :
    K2-99 :
    Kepler-1515 :
    Kepler-1517 :
    Kepler-1658 :
    Kepler-411 :
    Kepler-91 :
    KOI-12 :
    NGTS-24 :
    NGTS-9 :
    PH2 :
    TIC 237913194 :
    TOI-1107 :
    TOI-1181 :
    TOI-1194 :
    TOI-1246 :
    TOI-1338 : TOI-1338 A
    TOI-1408 :
    TOI-1416 :
    TOI-1420 :
    TOI-1452 :
    TOI-1516 :
    TOI-1680 :
    TOI-1811 :
    TOI-1820 :
    TOI-1937 A :
    TOI-198 :
    TOI-2025 :
    TOI-2046 :
    TOI-2048 :
    TOI-206 :
    TOI-2081 :
    TOI-2145 :
    TOI-2152 A : TOI-2152
    TOI-2154 :
    TOI-2158 :
    TOI-2193 A :
    TOI-2194 :
    TOI-2202 :
    TOI-2207 :
    TOI-2236 :
    TOI-2338 :
    TOI-2421 :
    TOI-2497 :
    TOI-2524 :
    TOI-2567 :
    TOI-2570 :
    TOI-2583 A :
    TOI-2587 A :
    TOI-262 :
    TOI-2641 :
    TOI-2796 :
    TOI-2803 A :
    TOI-2818 :
    TOI-2842 :
    TOI-2977 :
    TOI-3023 :
    TOI-3235 :
    TOI-3331 A :
    TOI-3364 :
    TOI-3540 A :
    TOI-3688 A :
    TOI-3693 :
    TOI-3714 :
    TOI-3807 :
    TOI-3819 :
    TOI-3884 :
    TOI-3912 :
    TOI-3976 A :
    TOI-4087 :
    TOI-411 : HD 22946
    TOI-4137 :
    TOI-4145 A :
    TOI-4308 :
    TOI-4342 :
    TOI-4377 :
    TOI-4406 :
    TOI-4463 A :
    TOI-4551 :
    TOI-4603 :
    TOI-4791 :
    TOI-615 :
    TOI-622 :
    TOI-778 :
    TOI-858 : TOI-858 B
    WASP-130 :
    WASP-150 :
    WASP-162 :
    WASP-171 :
    WASP-193 :
    WASP-60 :
    HD 110067 :
    TOI-3261 :
    CoRoT-10 :
    CoRoT-12 :
    CoRoT-13 :
    CoRoT-14 :
    CoRoT-16 :
    CoRoT-18 :
    CoRoT-22 :
    CoRoT-24 :
    CoRoT-25 :
    CoRoT-26 :
    CoRoT-28 :
    CoRoT-29 :
    CoRoT-31 :
    CoRoT-4 :
    CoRoT-6 :
    CoRoT-8 :
    CoRoT-9 :
    CoRoT-32 :
    EPIC 212297394 :
    EPIC 220674823 :
    EPIC 229004835 :
    EPIC 249893012 :
    HATS-14 :
    HATS-15 :
    HATS-16 :
    HATS-17 :
    HATS-32 :
    HATS-36 :
    HATS-44 :
    HATS-49 :
    HATS-54 :
    HATS-59 :
    HATS-63 :
    HATS-69 :
    HATS-71 :
    HATS-74 A :
    HATS-75 :
    HATS-76 :
    HATS-77 :
    HATS-8 :
    HD 137496 :
    HD 18599 :
    HD 20329 :
    HD 80653 :
    HD 88986 :
    K2-10 :
    K2-100 :
    K2-108 :
    K2-110 :
    K2-111 :
    K2-113 :
    K2-114 :
    K2-115 :
    K2-117 :
    K2-12 :
    K2-122 :
    K2-126 :
    K2-127 :
    K2-128 :
    K2-130 :
    K2-131 :
    K2-133 :
    K2-14 :
    K2-146 :
    K2-148 :
    K2-151 :
    K2-152 :
    K2-153 :
    K2-154 :
    K2-158 :
    K2-159 :
    K2-160 :
    K2-161 :
    K2-162 :
    K2-163 :
    K2-164 :
    K2-165 :
    K2-168 :
    K2-17 :
    K2-170 :
    K2-171 :
    K2-172 :
    K2-173 :
    K2-175 :
    K2-176 :
    K2-177 :
    K2-178 :
    K2-179 :
    K2-180 :
    K2-181 :
    K2-182 :
    K2-184 :
    K2-186 :
    K2-188 :
    K2-189 :
    K2-190 :
    K2-193 :
    K2-194 :
    K2-195 :
    K2-196 :
    K2-197 :
    K2-199 :
    K2-200 :
    K2-201 :
    K2-203 :
    K2-204 :
    K2-206 :
    K2-208 :
    K2-216 :
    K2-216 :
    K2-217 :
    K2-219 :
    K2-220 :
    K2-221 :
    K2-224 :
    K2-225 :
    K2-226 :
    K2-227 :
    K2-228 :
    K2-229 :
    K2-230 :
    K2-241 :
    K2-242 :
    K2-245 :
    K2-246 :
    K2-247 :
    K2-248 :
    K2-250 :
    K2-253 :
    K2-254 :
    K2-255 :
    K2-265 :
    K2-268 :
    K2-27 :
    K2-270 :
    K2-271 :
    K2-272 :
    K2-273 :
    K2-274 :
    K2-275 :
    K2-276 :
    K2-277 :
    K2-281 :
    K2-283 :
    K2-286 :
    K2-291 :
    K2-292 :
    K2-294 :
    K2-308 :
    K2-331 :
    K2-341 :
    K2-342 :
    K2-346 :
    K2-35 :
    K2-352 :
    K2-354 :
    K2-357 :
    K2-365 :
    K2-366 :
    K2-367 :
    K2-368 :
    K2-369 :
    K2-37 :
    K2-371 :
    K2-372 :
    K2-374 :
    K2-376 :
    K2-378 :
    K2-379 :
    K2-38 :
    K2-380 :
    K2-381 :
    K2-383 :
    K2-387 :
    K2-388 :
    K2-389 :
    K2-390 :
    K2-391 :
    K2-393 :
    K2-394 :
    K2-395 :
    K2-396 :
    K2-397 :
    K2-398 :
    K2-401 :
    K2-404 :
    K2-407 :
    K2-408 :
    K2-409 :
    K2-416 :
    K2-43 :
    K2-44 :
    K2-45 :
    K2-5 :
    K2-62 :
    K2-65 :
    K2-68 :
    K2-70 :
    K2-73 :
    K2-75 :
    K2-77 :
    K2-8 :
    K2-80 :
    K2-81 :
    K2-84 :
    K2-86 :
    K2-98 :
    KIC 9663113 :
    KOI-1257 :
    KOI-217 :
    KOI-351 :
    KOI-3680 :
    KOI-7368 :
    KOI-984 :
    Kepler-1002 :
    Kepler-1003 :
    Kepler-1004 :
    Kepler-101 :
    Kepler-1011 :
    Kepler-1012 :
    Kepler-1016 :
    Kepler-1019 :
    Kepler-1028 :
    Kepler-103 :
    Kepler-1035 :
    Kepler-1037 :
    Kepler-1039 :
    Kepler-1055 :
    Kepler-1061 :
    Kepler-1065 :
    Kepler-1066 :
    Kepler-107 :
    Kepler-1072 :
    Kepler-1073 :
    Kepler-1075 :
    Kepler-1082 :
    Kepler-109 :
    Kepler-1091 :
    Kepler-1099 :
    Kepler-1106 :
    Kepler-1107 :
    Kepler-1139 :
    Kepler-1148 :
    Kepler-1167 :
    Kepler-117 :
    Kepler-1179 :
    Kepler-1181 :
    Kepler-1189 :
    Kepler-1193 :
    Kepler-1205 :
    Kepler-1206 :
    Kepler-1217 :
    Kepler-1228 :
    Kepler-1256 :
    Kepler-1259 :
    Kepler-1270 :
    Kepler-1271 :
    Kepler-1310 :
    Kepler-1311 :
    Kepler-1312 :
    Kepler-1313 :
    Kepler-1315 :
    Kepler-1317 :
    Kepler-1320 :
    Kepler-1323 :
    Kepler-1344 :
    Kepler-1346 :
    Kepler-1356 :
    Kepler-148 :
    Kepler-15 :
    Kepler-1514 :
    Kepler-1519 :
    Kepler-1520 :
    Kepler-1522 :
    Kepler-1523 :
    Kepler-1524 :
    Kepler-1528 :
    Kepler-1530 :
    Kepler-154 :
    Kepler-157 :
    Kepler-1604 :
    Kepler-1613 :
    Kepler-1624 :
    Kepler-1627 :
    Kepler-1642 :
    Kepler-1654 :
    Kepler-1655 :
    Kepler-166 :
    Kepler-1704 :
    Kepler-198 :
    Kepler-289 :
    Kepler-30 :
    Kepler-323 :
    Kepler-4 :
    Kepler-41 :
    Kepler-412 :
    Kepler-418 :
    Kepler-421 :
    Kepler-422 :
    Kepler-423 :
    Kepler-424 :
    Kepler-425 :
    Kepler-426 :
    Kepler-427 :
    Kepler-428 :
    Kepler-43 :
    Kepler-44 :
    Kepler-46 :
    Kepler-460 : KIC 5437945
    Kepler-463 :
    Kepler-464 :
    Kepler-466 :
    Kepler-472 :
    Kepler-473 :
    Kepler-474 :
    Kepler-475 :
    Kepler-476 :
    Kepler-485 :
    Kepler-487 :
    Kepler-490 :
    Kepler-491 :
    Kepler-495 :
    Kepler-497 :
    Kepler-498 :
    Kepler-501 :
    Kepler-502 :
    Kepler-507 :
    Kepler-511 :
    Kepler-514 :
    Kepler-518 :
    Kepler-523 :
    Kepler-529 :
    Kepler-530 :
    Kepler-531 :
    Kepler-532 :
    Kepler-533 :
    Kepler-535 :
    Kepler-536 :
    Kepler-537 :
    Kepler-539 :
    Kepler-541 :
    Kepler-543 :
    Kepler-546 :
    Kepler-547 :
    Kepler-548 :
    Kepler-550 :
    Kepler-552 :
    Kepler-553 :
    Kepler-554 :
    Kepler-557 :
    Kepler-56 :
    Kepler-561 :
    Kepler-562 :
    Kepler-564 :
    Kepler-565 :
    Kepler-568 :
    Kepler-578 :
    Kepler-585 :
    Kepler-586 :
    Kepler-590 :
    Kepler-592 :
    Kepler-605 :
    Kepler-608 :
    Kepler-609 :
    Kepler-611 :
    Kepler-612 :
    Kepler-618 :
    Kepler-619 :
    Kepler-620 :
    Kepler-621 :
    Kepler-629 :
    Kepler-634 :
    Kepler-636 :
    Kepler-645 :
    Kepler-648 :
    Kepler-650 :
    Kepler-656 :
    Kepler-664 :
    Kepler-666 :
    Kepler-667 :
    Kepler-669 :
    Kepler-670 :
    Kepler-673 :
    Kepler-674 :
    Kepler-677 :
    Kepler-678 :
    Kepler-680 :
    Kepler-682 :
    Kepler-683 :
    Kepler-684 :
    Kepler-685 :
    Kepler-686 :
    Kepler-688 :
    Kepler-690 :
    Kepler-693 :
    Kepler-694 :
    Kepler-695 :
    Kepler-696 :
    Kepler-697 :
    Kepler-700 :
    Kepler-702 :
    Kepler-703 :
    Kepler-707 :
    Kepler-708 :
    Kepler-714 :
    Kepler-715 :
    Kepler-718 :
    Kepler-719 :
    Kepler-720 :
    Kepler-722 :
    Kepler-723 :
    Kepler-724 :
    Kepler-725 :
    Kepler-728 :
    Kepler-729 :
    Kepler-731 :
    Kepler-734 :
    Kepler-736 :
    Kepler-74 :
    Kepler-740 :
    Kepler-743 :
    Kepler-75 :
    Kepler-750 :
    Kepler-755 :
    Kepler-762 :
    Kepler-767 :
    Kepler-77 :
    Kepler-773 :
    Kepler-775 :
    Kepler-78 :
    Kepler-782 :
    Kepler-785 :
    Kepler-796 :
    Kepler-799 :
    Kepler-80 :
    Kepler-805 :
    Kepler-808 :
    Kepler-814 :
    Kepler-815 :
    Kepler-816 :
    Kepler-817 :
    Kepler-818 :
    Kepler-822 :
    Kepler-825 :
    Kepler-827 :
    Kepler-828 :
    Kepler-831 :
    Kepler-842 :
    Kepler-843 :
    Kepler-845 :
    Kepler-849 :
    Kepler-855 :
    Kepler-856 :
    Kepler-857 :
    Kepler-858 :
    Kepler-860 :
    Kepler-867 :
    Kepler-872 :
    Kepler-890 :
    Kepler-891 :
    Kepler-893 :
    Kepler-904 :
    Kepler-905 :
    Kepler-908 :
    Kepler-912 :
    Kepler-914 :
    Kepler-915 :
    Kepler-922 :
    Kepler-932 :
    Kepler-943 :
    Kepler-950 :
    Kepler-951 :
    Kepler-952 :
    Kepler-953 :
    Kepler-954 :
    Kepler-956 :
    Kepler-957 :
    Kepler-960 :
    Kepler-972 :
    Kepler-975 :
    Kepler-990 :
    Kepler-994 :
    Kepler-996 :
    Kepler-997 :
    LHS 1815 :
    NGTS-1 :
    NGTS-15 :
    NGTS-16 :
    NGTS-17 :
    NGTS-18 :
    NGTS-20 :
    NGTS-21 :
    NGTS-23 :
    NGTS-25 :
    NGTS-3 A :
    NGTS-4 :
    OGLE-TR-182 :
    POTS-1 :
    TIC 172900988 : TIC 172900988 Aa
    TIC 279401253 :
    TOI-1052 :
    TOI-1062 :
    TOI-1221 :
    TOI-1235 :
    TOI-1444 :
    TOI-1696 :
    TOI-1736 :
    TOI-2084 :
    TOI-2095 :
    TOI-2096 :
    TOI-2141 :
    TOI-2180 :
    TOI-2184 :
    TOI-2196 :
    TOI-2257 :
    TOI-2285 :
    TOI-244 :
    TOI-2525 :
    TOI-2589 :
    TOI-3757 :
    TOI-3984 A :
    TOI-4127 :
    TOI-4184 :
    TOI-4201 :
    TOI-4562 :
    TOI-4582 :
    TOI-4860 :
    TOI-5205 :
    TOI-5293 A :
    TOI-5542 :
    TOI-733 :
    TOI-763 :
    TOI-784 : HD 307842
    WASP-59 :
    WTS-2 :
    Wendelstein-1 :
    Wendelstein-2 :
    TOI-6086 :
    GJ 12 : Gliese 12
    TOI-1135 :
    TOI-4336 A :
    SPECULOOS-3 :
    TOI-904 :
    TOI-4559 :
    TOI-771 :
    TOI-6008 :
    TOI-5799 :
    TOI-2120 :
    TOI-782 :
    TOI-260 :
    BD+05 4868 : BD+05 4868 A
    HD 12572 : HIP 9618
    LP 890-9 :
    TOI-6255 :
    TOI-6894 :
    GJ 341 :
    BD+20 594 :
    BD-14 3065 A :
    CoRoT-20 :
    CoRoT-32 :
    EPIC 201595106 :
    EPIC 206024342 :
    EPIC 212587672 :
    EPIC 212624936 :
    EPIC 212737443 :
    K2-282 :
    HAT-P-63 :
    HATS-66 :
    HIP 56998 :
    HD 114082 :
    HD 135694 :
    HD 158259 :
    HD 21520 :
    HD 224018 :
    HD 25463 :
    HD 35843 :
    HD 56414 :
    HD 6061 :
    HD 73344 :
    HD 77946 :
    HIP 8152 :
    IRAS 04125+2902 :
    K2-101 :
    K2-118 :
    K2-119 :
    K2-123 :
    K2-13 :
    K2-137 :
    K2-147 :
    K2-149 :
    K2-15 :
    K2-150 :
    K2-156 :
    K2-16 :
    K2-183 :
    K2-185 :
    K2-187 :
    K2-202 :
    K2-205 :
    K2-207 :
    K2-209 :
    K2-214 :
    K2-215 :
    K2-218 :
    K2-231 :
    K2-243 :
    K2-244 :
    K2-249 :
    K2-251 :
    K2-252 :
    K2-258 :
    K2-263 :
    K2-278 :
    K2-279 :
    K2-290 :
    K2-318 :
    K2-319 :
    K2-322 :
    K2-323 :
    K2-324 :
    K2-330 :
    K2-337 :
    K2-345 :
    K2-349 :
    K2-350 :
    K2-355 :
    K2-356 :
    K2-358 :
    K2-373 :
    K2-382 :
    K2-384 :
    K2-385 :
    K2-386 :
    K2-4 :
    K2-400 :
    K2-402 :
    K2-414 :
    K2-42 :
    K2-48 :
    K2-50 :
    K2-54 :
    K2-57 :
    K2-59 :
    K2-6 :
    K2-61 :
    K2-63 :
    K2-64 :
    K2-69 :
    K2-7 :
    K2-72 :
    K2-74 :
    K2-83 :
    K2-85 :
    K2-88 :
    K2-91 :
    KIC 3558849 :
    KOI-7892 :
    KOI-134 :
    KOI-142 :
    KOI-1783 :
    Kepler-968 :
    KOI-7913 A :
    Kepler-100 :
    Kepler-106 :
    Kepler-112 :
    Kepler-113 :
    Kepler-114 :
    Kepler-116 :
    Kepler-1171 :
    Kepler-119 :
    Kepler-122 :
    Kepler-130 :
    Kepler-131 :
    Kepler-1311 :
    Kepler-1312 :
    Kepler-1313 :
    Kepler-1326 :
    Kepler-134 :
    Kepler-138 :
    Kepler-139 :
    Kepler-142 :
    Kepler-1442 :
    Kepler-146 :
    Kepler-149 :
    Kepler-1518 :
    Kepler-152 :
    Kepler-1521 :
    Kepler-153 :
    Kepler-155 :
    Kepler-161 :
    Kepler-162 :
    Kepler-1643 :
    Kepler-1661 :
    Kepler-1665 :
    Kepler-1672 :
    Kepler-1676 :
    Kepler-170 :
    Kepler-1709 :
    Kepler-1710 :
    Kepler-1712 :
    Kepler-1714 :
    Kepler-1715 :
    Kepler-1716 :
    Kepler-1718 :
    Kepler-1719 :
    Kepler-1725 :
    Kepler-1726 :
    Kepler-1731 :
    Kepler-1746 :
    Kepler-1755 :
    Kepler-1769 :
    Kepler-1771 :
    Kepler-1772 :
    Kepler-1782 :
    Kepler-1791 :
    Kepler-1832 :
    Kepler-186 :
    Kepler-1868 :
    Kepler-1869 :
    Kepler-1870 :
    Kepler-192 :
    Kepler-1929 :
    Kepler-1935 :
    Kepler-1939 :
    Kepler-199 :
    Kepler-1990 :
    Kepler-202 :
    Kepler-203 :
    Kepler-207 :
    Kepler-21 :
    Kepler-210 :
    Kepler-220 :
    Kepler-221 :
    Kepler-261 :
    Kepler-310 :
    Kepler-314 :
    Kepler-318 :
    Kepler-319 :
    Kepler-324 :
    Kepler-350 :
    Kepler-36 :
    Kepler-410 A :
    Kepler-414 :
    Kepler-433 :
    Kepler-449 :
    Kepler-453 :
    Kepler-461 :
    Kepler-462 :
    Kepler-465 :
    Kepler-471 :
    Kepler-477 :
    Kepler-478 :
    Kepler-479 :
    Kepler-480 :
    Kepler-484 :
    Kepler-499 :
    Kepler-50 :
    Kepler-506 :
    Kepler-509 :
    Kepler-510 :
    Kepler-516 :
    Kepler-517 :
    Kepler-519 :
    Kepler-522 :
    Kepler-538 :
    Kepler-560 :
    Kepler-569 :
    Kepler-572 :
    Kepler-617 :
    Kepler-622 :
    Kepler-643 :
    Kepler-65 :
    Kepler-652 :
    Kepler-69 :
    Kepler-732 :
    Kepler-753 :
    Kepler-783 :
    Kepler-8 :
    Kepler-803 :
    Kepler-820 :
    Kepler-880 :
    Kepler-901 :
    Kepler-949 :
    Kepler-95 :
    Kepler-959 :
    Kepler-971 :
    Kepler-974 :
    Kepler-98 :
    NGTS-27 :
    NGTS-30 :
    NGTS-31 :
    NGTS-32 :
    NGTS-33 :
    Ross 176 :
    TIC 139270665 :
    TIC 365102760 :
    TIC 393818343 :
    TIC 434398831 :
    TIC 46432937 :
    TIC 88785435 :
    TOI-1011 :
    TOI-1117 :
    TOI-1173 :
    TOI-1174 :
    TOI-1180 :
    TOI-1184 :
    TOI-1199 :
    TOI-1203 :
    TOI-1224 :
    TOI-1238 :
    TOI-1244 :
    TOI-1248 :
    TOI-1249 :
    TOI-1269 :
    TOI-1273 :
    TOI-1279 :
    TOI-128 :
    TOI-1294 :
    TOI-1295 :
    TOI-1301 :
    TOI-1346 :
    TOI-1347 :
    TOI-1386 :
    TOI-1410 :
    TOI-1437 :
    TOI-1438 :
    TOI-1439 :
    TOI-1443 :
    TOI-1448 :
    TOI-1450 A :
    TOI-1451 :
    TOI-1453 :
    TOI-1467 :
    TOI-1472 :
    TOI-1630 :
    TOI-1659 :
    TOI-1669 :
    TOI-1683 :
    TOI-1691 :
    TOI-1716 :
    TOI-1723 :
    TOI-1739 :
    TOI-1742 :
    TOI-1743 :
    TOI-1744 :
    TOI-1749 :
    TOI-1751 :
    TOI-1753 :
    TOI-1758 :
    TOI-1768 :
    TOI-1772 :
    TOI-1775 :
    TOI-1776 :
    TOI-1777 :
    TOI-1782 :
    TOI-1794 :
    TOI-1798 :
    TOI-1799 :
    TOI-1803 :
    TOI-1806 :
    TOI-1823 :
    TOI-1824 :
    TOI-1836 :
    TOI-1846 :
    TOI-1855 :
    TOI-1883 :
    TOI-1898 :
    TOI-2005 :
    TOI-2015 :
    TOI-2019 :
    TOI-2031 A :
    TOI-2068 :
    TOI-2088 :
    TOI-2107 :
    TOI-2128 :
    TOI-2169 A :
    TOI-2211 :
    TOI-2266 :
    TOI-2274 :
    TOI-2295 :
    TOI-2322 :
    TOI-2328 :
    TOI-2346 :
    TOI-2368 :
    TOI-2374 :
    TOI-2379 :
    TOI-238 :
    TOI-2382 :
    TOI-2384 :
    TOI-2407 :
    TOI-2420 :
    TOI-2447 :
    TOI-2449 :
    TOI-2458 :
    TOI-2485 :
    TOI-2529 :
    TOI-2537 :
    TOI-2545 :
    TOI-2580 :
    TOI-260 :
    TOI-261 :
    TOI-2714 :
    TOI-2719 :
    TOI-2768 :
    TOI-286 :
    TOI-2876 :
    TOI-2886 :
    TOI-2969 :
    TOI-2981 :
    TOI-2986 :
    TOI-2989 :
    TOI-2992 :
    TOI-3071 :
    TOI-3135 :
    TOI-3160 A :
    TOI-329 :
    TOI-3321 :
    TOI-3353 :
    TOI-3464 :
    TOI-3474 :
    TOI-3486 :
    TOI-3493 :
    TOI-3523 A :
    TOI-3568 :
    TOI-3593 :
    TOI-3682 :
    TOI-3837 :
    TOI-3856 :
    TOI-3877 :
    TOI-3894 :
    TOI-3919 :
    TOI-3980 :
    TOI-406 :
    TOI-4153 :
    TOI-4155 :
    TOI-4214 :
    TOI-4364 :
    TOI-4379 :
    TOI-4438 :
    TOI-4443 :
    TOI-4465 :
    TOI-4487 A :
    TOI-4495 :
    TOI-4504 :
    TOI-4515 :
    TOI-4527 :
    TOI-4602 :
    TOI-4638 :
    TOI-4734 :
    TOI-4773 :
    TOI-4794 :
    TOI-480 :
    TOI-4914 :
    TOI-4961 :
    TOI-4994 :
    TOI-5005 :
    TOI-5027 :
    TOI-5076 :
    TOI-5082 :
    TOI-5108 :
    TOI-5110 :
    TOI-512 :
    TOI-5143 :
    TOI-5153 :
    TOI-5174 :
    TOI-5181 A :
    TOI-5210 :
    TOI-5232 :
    TOI-5238 :
    TOI-5300 :
    TOI-5301 :
    TOI-5319 :
    TOI-5322 :
    TOI-5340 :
    TOI-5350 :
    TOI-5386 A :
    TOI-5388 :
    TOI-5573 :
    TOI-558 :
    TOI-5592 :
    TOI-5713 :
    TOI-5720 :
    TOI-5726 :
    TOI-5786 :
    TOI-5795 :
    TOI-5800 :
    TOI-5817 :
    TOI-6000 :
    TOI-6002 :
    TOI-6016 :
    TOI-6029 :
    TOI-6034 :
    TOI-6038 A :
    TOI-6054 :
    TOI-6109 :
    TOI-6130 :
    TOI-6223 :
    TOI-6303 :
    TOI-6324 :
    TOI-6330 :
    TOI-6478 :
    TOI-654 :
    TOI-6628 :
    TOI-663 :
    TOI-6651 :
    TOI-669 :
    TOI-6695 :
    TOI-7041 :
    TOI-715 :
    TOI-757 :
    TOI-762 A :
    TOI-815 :
    TOI-871 :
    TOI-880 :
    TOI-907 :
    WASP-102 :
    WASP-116 :
    WASP-149 :
    WASP-154 :
    WASP-155 :
    WASP-188 :
    WASP-194 :
    WASP-195 :
    WASP-197 :
    testJup :
    '''
    # BD+05 4868   BD+05 4868 A
    # HD 12572     HIP 9618
    # LP-890-9     rename LP 890-9
    # TOI-2031        TOI-2031.01 candidate, so not in main table
    # TOI-2431        TOI-2431.01 candidate, so not in main table
    # TOI-6255     ok
    # TOI-6894     ok
    # TYC-7052-1753-1 TOI-2490.01 candidate, so not in main table

    # these two JWST targets are not yet listed in the Exoplanet Archive composite table:
    #  6/21/25 first one is ok now!  (2025 paper)  no HIP name though
    # GJ 341 : HIP 45908  (12/16/24 still listed as a candidate planet)
    # GJ 1008 : HIP 1532  (12/16/24 OK this one is there now, as TOI-260)

    return
