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
    CoRoTID 223977153 :
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
    testJup :
    BD+05 4868 : BD+05 4868 A
    HD 12572 : HIP 9618
    LP 890-9 :
    TOI-6255 :
    TOI-6894 :
    GJ 341 :
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


# ----------------- --------------------------------------------------
# -- TARGET ON DISK -- -----------------------------------------------
# FIRST COL MATCHES TARGET LIST NAME
# SECOND COL IS THE DATA FOLDER IN /proj/sdp/data/sci
def targetondisk():
    '''
    55 Cnc : 55CNC
    CoRoT-1 : COROT1
    CoRoT-2 : COROT2
    CoRoT-7 : COROT7
    GJ 1132 : GJ1132
    GJ 1214 : GJ1214
    GJ 1252 : GJ1252
    GJ 3053 : GJ3053
    GJ 3470 : GJ3470
    GJ 436 : GJ436
    GJ 9827 : GJ9827
    HAT-P-1 : HATP1
    HAT-P-11 : HATP11
    HAT-P-12 : HATP12
    HAT-P-13 : HATP13
    HAT-P-15 : HATP15
    HAT-P-16 : HATP16
    HAT-P-17 : HATP17
    HAT-P-18 : HATP18
    HAT-P-19 : HATP19
    HAT-P-2 : HATP2
    HAT-P-20 : HATP20
    HAT-P-22 : HATP22
    HAT-P-23 : HATP23
    HAT-P-26 : HATP26
    HAT-P-3 : HATP3
    HAT-P-30 : HATP30
    HAT-P-32 : HATP32
    HAT-P-33 : HATP33
    HAT-P-34 : HATP34
    HAT-P-38 : HATP38
    HAT-P-4 : HATP4
    HAT-P-40 : HATP40
    HAT-P-41 : HATP41
    HAT-P-5 : HATP5
    HAT-P-6 : HATP6
    HAT-P-7 : HATP7
    HAT-P-8 : HATP8
    HATS-28 : HATS28
    HATS-3 : HATS3
    HATS-65 : HATS65
    HATS-7 : HATS7
    HD 106315 : HD106315
    HD 149026 : HD149026
    HD 17156 : HD17156
    HD 189733 : HD189733
    HD 209458 : HD209458
    HD 213885 : HD213885
    HD 219134 : HD219134
    HD 219666 : HD219666
    HD 23472 : HD23472
    HD 97658 : HD97658
    HR 858 : HR858
    K2-124 : K2124
    K2-132 : K2132
    K2-136 : K2136
    K2-138 : K2138
    K2-141 : K2141
    K2-167 : K2167
    K2-174 : K2174
    K2-18 : K218
    K2-19 : K219
    K2-21 : K221
    K2-212 : K2212
    K2-22 : K222
    K2-24 : K224
    K2-25 : K225
    K2-26 : K226
    K2-28 : K228
    K2-289 : K2289
    K2-3 : K23
    K2-31 : K231
    K2-32 : K232
    K2-33 : K233
    K2-36 : K236
    K2-52 : K252
    K2-53 : K253
    K2-55 : K255
    K2-58 : K258
    K2-79 : K279
    K2-87 : K287
    K2-9 : K29
    K2-90 : K290
    K2-93 : K293
    K2-95 : K295
    K2-96 : K296
    K2-97 : K297
    KELT-1 : KELT1
    KELT-11 : KELT11
    KELT-14 : KELT14
    KELT-16 : KELT16
    KELT-20 : KELT20
    KELT-3 : KELT3
    KELT-7 : KELT7
    KELT-9 : KELT9
    Kepler-10 : KEPLER10
    Kepler-11 : KEPLER11
    Kepler-102 : KEPLER102
    Kepler-104 : KEPLER104
    Kepler-1083 : KEPLER1083
    Kepler-12 : KEPLER12
    Kepler-125 : KEPLER125
    Kepler-126 : KEPLER126
    Kepler-127 : KEPLER127
    Kepler-13 : KEPLER13
    Kepler-1339 : KEPLER1339
    Kepler-138 : KEPLER138
    Kepler-14 : KEPLER14
    Kepler-1485 : KEPLER1485
    Kepler-1492 : KEPLER1492
    Kepler-156 : KEPLER156
    Kepler-1568 : KEPLER1568
    Kepler-158 : KEPLER158
    Kepler-16 : KEPLER16
    Kepler-1625 : KEPLER1625
    Kepler-1651 : KEPLER1651
    Kepler-167 : KEPLER167
    Kepler-17 : KEPLER17
    Kepler-18 : KEPLER18
    Kepler-19 : KEPLER19
    Kepler-20 : KEPLER20
    Kepler-205 : KEPLER205
    Kepler-218 : KEPLER218
    Kepler-236 : KEPLER236
    Kepler-249 : KEPLER249
    Kepler-25 : KEPLER25
    Kepler-26 : KEPLER26
    Kepler-293 : KEPLER293
    Kepler-297 : KEPLER297
    Kepler-309 : KEPLER309
    Kepler-32 : KEPLER32
    Kepler-37 : KEPLER37
    Kepler-395 : KEPLER395
    Kepler-45 : KEPLER45
    Kepler-454 : KEPLER454
    Kepler-48 : KEPLER48
    Kepler-482 : KEPLER482
    Kepler-49 : KEPLER49
    Kepler-5 : KEPLER5
    Kepler-504 : KEPLER504
    Kepler-505 : KEPLER505
    Kepler-570 : KEPLER570
    Kepler-582 : KEPLER582
    Kepler-598 : KEPLER598
    Kepler-6 : KEPLER6
    Kepler-603 : KEPLER603
    Kepler-61 : KEPLER61
    Kepler-62 : KEPLER62
    Kepler-68 : KEPLER68
    Kepler-7 : KEPLER7
    Kepler-705 : KEPLER705
    Kepler-737 : KEPLER737
    Kepler-769 : KEPLER769
    Kepler-786 : KEPLER786
    Kepler-9 : KEPLER9
    Kepler-93 : KEPLER93
    Kepler-94 : KEPLER94
    LHS 3844 : LHS3844
    OGLE-TR-056 : OGLETR056
    OGLE-TR-10 : OGLETR10
    Qatar-1 : QATAR1
    Qatar-2 : QATAR2
    TOI-270 : TOI270
    TOI-700 : TOI700
    TOI-849 : TOI849
    TRAPPIST-1 : TRAPPIST1
    TrES-1 : TRES1
    TrES-2 : TRES2
    TrES-3 : TRES3
    WASP-1 : WASP1
    WASP-10 : WASP10
    WASP-100 : WASP100
    WASP-101 : WASP101
    WASP-103 : WASP103
    WASP-104 : WASP104
    WASP-107 : WASP107
    WASP-11 : WASP11
    WASP-12 : WASP12
    WASP-121 : WASP121
    WASP-127 : WASP127
    WASP-13 : WASP13
    WASP-131 : WASP131
    WASP-14 : WASP14
    WASP-140 : WASP140
    WASP-15 : WASP15
    WASP-16 : WASP16
    WASP-17 : WASP17
    WASP-18 : WASP18
    WASP-19 : WASP19
    WASP-2 : WASP2
    WASP-21 : WASP21
    WASP-24 : WASP24
    WASP-26 : WASP26
    WASP-28 : WASP28
    WASP-29 : WASP29
    WASP-3 : WASP3
    WASP-31 : WASP31
    WASP-32 : WASP32
    WASP-33 : WASP33
    WASP-34 : WASP34
    WASP-35 : WASP35
    WASP-36 : WASP36
    WASP-37 : WASP37
    WASP-38 : WASP38
    WASP-39 : WASP39
    WASP-4 : WASP4
    WASP-43 : WASP43
    WASP-46 : WASP46
    WASP-48 : WASP48
    WASP-49 : WASP49
    WASP-5 : WASP5
    WASP-50 : WASP50
    WASP-52 : WASP52
    WASP-6 : WASP6
    WASP-62 : WASP62
    WASP-63 : WASP63
    WASP-64 : WASP64
    WASP-65 : WASP65
    WASP-67 : WASP67
    WASP-69 : WASP69
    WASP-7 : WASP7
    WASP-72 : WASP72
    WASP-74 : WASP74
    WASP-75 : WASP75
    WASP-76 : WASP76
    WASP-77 : WASP77A
    WASP-78 : WASP78
    WASP-79 : WASP79
    WASP-8 : WASP8
    WASP-80 : WASP80
    WASP-87 : WASP87
    WASP-94 : WASP94A
    WASP-95 : WASP95
    WASP-96 : WASP96
    WASP-97 : WASP97
    XO-1 : XO1
    XO-2 : XO2
    XO-3 : XO3
    XO-4 : XO4
    XO-5 : XO5
    AU Mic : AUMIC
    CoRoT-5 : COROT5
    CoRoT-11 : COROT11
    CoRoT-19 : COROT19
    DS Tuc A : DSTUCA
    EPIC 211945201 : EPIC211945201
    EPIC 246851721 : EPIC246851721
    G 9-40 : G940
    GJ 3473 : GJ3473
    GJ 357 : GJ357
    GJ 367 : GJ367
    GJ 3929 : GJ3929
    GJ 486 : GJ486
    GPX-1 : GPX1
    HAT-P-14 : HATP14
    HAT-P-21 : HATP21
    HAT-P-24 : HATP24
    HAT-P-25 : HATP25
    HAT-P-27 : HATP27
    HAT-P-28 : HATP28
    HAT-P-29 : HATP29
    HAT-P-31 : HATP31
    HAT-P-35 : HATP35
    HAT-P-36 : HATP36
    HAT-P-37 : HATP37
    HAT-P-39 : HATP39
    HAT-P-42 : HATP42
    HAT-P-43 : HATP43
    HAT-P-44 : HATP44
    HAT-P-45 : HATP45
    HAT-P-46 : HATP46
    HAT-P-49 : HATP49
    HAT-P-50 : HATP50
    HAT-P-51 : HATP51
    HAT-P-52 : HATP52
    HAT-P-53 : HATP53
    HAT-P-54 : HATP54
    HAT-P-55 : HATP55
    HAT-P-56 : HATP56
    HAT-P-57 : HATP57
    HAT-P-58 : HATP58
    HAT-P-59 : HATP59
    HAT-P-60 : HATP60
    HAT-P-61 : HATP61
    HAT-P-62 : HATP62
    HAT-P-64 : HATP64
    HAT-P-65 : HATP65
    HAT-P-66 : HATP66
    HAT-P-67 : HATP67
    HAT-P-68 : HATP68
    HAT-P-69 : HATP69
    HAT-P-70 : HATP70
    HAT-P-9 : HATP9
    HATS-1 : HATS1
    HATS-11 : HATS11
    HATS-13 : HATS13
    HATS-18 : HATS18
    HATS-2 : HATS2
    HATS-23 : HATS23
    HATS-24 : HATS24
    HATS-25 : HATS25
    HATS-26 : HATS26
    HATS-27 : HATS27
    HATS-29 : HATS29
    HATS-30 : HATS30
    HATS-31 : HATS31
    HATS-33 : HATS33
    HATS-34 : HATS34
    HATS-35 : HATS35
    HATS-37 A : HATS37A
    HATS-38 : HATS38
    HATS-39 : HATS39
    HATS-4 : HATS4
    HATS-40 : HATS40
    HATS-41 : HATS41
    HATS-42 : HATS42
    HATS-43 : HATS43
    HATS-46 : HATS46
    HATS-47 : HATS47
    HATS-48 A : HATS48A
    HATS-5 : HATS5
    HATS-50 : HATS50
    HATS-51 : HATS51
    HATS-52 : HATS52
    HATS-53 : HATS53
    HATS-56 : HATS56
    HATS-57 : HATS57
    HATS-58 A : HATS58A
    HATS-6 : HATS6
    HATS-60 : HATS60
    HATS-62 : HATS62
    HATS-64 : HATS64
    HATS-67 : HATS67
    HATS-68 : HATS68
    HATS-70 : HATS70
    HATS-72 : HATS72
    HATS-9 : HATS9
    HD 108236 : HD108236
    HD 110082 : HD110082
    HD 110113 : HD110113
    HD 136352 : HD136352
    HD 1397 : HD1397
    HD 152843 : HD152843
    HD 15337 : HD15337
    HD 183579 : HD183579
    HD 191939 : HD191939
    HD 202772 A : HD202772A
    HD 207897 : HD207897
    HD 221416 : HD221416
    HD 2685 : HD2685
    HD 332231 : HD332231
    HD 5278 : HD5278
    HD 63433 : HD63433
    HD 63935 : HD63935
    HD 73583 : HD73583
    HD 86226 : HD86226
    HD 89345 : HD89345
    HIP 65 A : HIP65A
    HIP 67522 : HIP67522
    K2-107 : K2107
    K2-116 : K2116
    K2-121 : K2121
    K2-129 : K2129
    K2-139 : K2139
    K2-140 : K2140
    K2-155 : K2155
    K2-198 : K2198
    K2-222 : K2222
    K2-232 : K2232
    K2-237 : K2237
    K2-238 : K2238
    K2-239 : K2239
    K2-260 : K2260
    K2-261 : K2261
    K2-266 : K2266
    K2-280 : K2280
    K2-284 : K2284
    K2-287 : K2287
    K2-29 : K229
    K2-295 : K2295
    K2-329 : K2329
    K2-333 : K2333
    K2-334 : K2334
    K2-34 : K234
    K2-353 : K2353
    K2-39 : K239
    K2-403 : K2403
    K2-405 : K2405
    K2-406 : K2406
    KELT-10 : KELT10
    KELT-12 : KELT12
    KELT-15 : KELT15
    KELT-17 : KELT17
    KELT-18 : KELT18
    KELT-19 A : KELT19A
    KELT-2 A : KELT2A
    KELT-21 : KELT21
    KELT-23 A : KELT23A
    KELT-24 : KELT24
    KELT-4 A : KELT4A
    KELT-6 : KELT6
    KELT-8 : KELT8
    KOI-13 : KOI13
    KOI-94 : KOI94
    KPS-1 : KPS1
    Kepler-105 : KEPLER105
    Kepler-108 : KEPLER108
    Kepler-1314 : KEPLER1314
    Kepler-1513 : KEPLER1513
    Kepler-33 : KEPLER33
    Kepler-396 : KEPLER396
    Kepler-42 : KEPLER42
    Kepler-435 : KEPLER435
    Kepler-444 : KEPLER444
    Kepler-447 : KEPLER447
    Kepler-450 : KEPLER450
    Kepler-468 : KEPLER468
    Kepler-489 : KEPLER489
    Kepler-76 : KEPLER76
    Kepler-79 : KEPLER79
    L 98-59 : L9859
    LHS 1478 : LHS1478
    LHS 1678 : LHS1678
    LP 714-47 : LP71447
    LP 791-18 : LP79118
    LTT 1445 A : LTT1445A
    LTT 3780 : LTT3780
    LTT 9779 : LTT9779
    MASCARA-1 : MASCARA1
    MASCARA-4 : MASCARA4
    NGTS-10 : NGTS10
    NGTS-11 : NGTS11
    NGTS-12 : NGTS12
    NGTS-13 : NGTS13
    NGTS-2 : NGTS2
    NGTS-5 : NGTS5
    NGTS-6 : NGTS6
    NGTS-8 : NGTS8
    Qatar-10 : QATAR10
    Qatar-4 : QATAR4
    Qatar-5 : QATAR5
    Qatar-6 : QATAR6
    Qatar-7 : QATAR7
    Qatar-8 : QATAR8
    Qatar-9 : QATAR9
    TIC 257060897 : TIC257060897
    TOI-1064 : TOI1064
    TOI-1075 : TOI1075
    TOI-1130 : TOI1130
    TOI-1201 : TOI1201
    TOI-122 : TOI122
    TOI-1227 : TOI1227
    TOI-1231 : TOI1231
    TOI-125 : TOI125
    TOI-1259 A : TOI1259A
    TOI-1260 : TOI1260
    TOI-1266 : TOI1266
    TOI-1268 : TOI1268
    TOI-1296 : TOI1296
    TOI-1298 : TOI1298
    TOI-1333 : TOI1333
    TOI-1411 : TOI1411
    TOI-1431 : TOI1431
    TOI-1442 : TOI1442
    TOI-1478 : TOI1478
    TOI-150 : TOI150
    TOI-1518 : TOI1518
    TOI-157 : TOI157
    TOI-1601 : TOI1601
    TOI-163 : TOI163
    TOI-1670 : TOI1670
    TOI-1685 : TOI1685
    TOI-169 : TOI169
    TOI-1693 : TOI1693
    TOI-172 : TOI172
    TOI-1728 : TOI1728
    TOI-1759 : TOI1759
    TOI-178 : TOI178
    TOI-1789 : TOI1789
    TOI-1807 : TOI1807
    TOI-1842 : TOI1842
    TOI-1860 : TOI1860
    TOI-1899 : TOI1899
    TOI-201 : TOI201
    TOI-2076 : TOI2076
    TOI-2109 : TOI2109
    TOI-216 : TOI216
    TOI-2260 : TOI2260
    TOI-2337 : TOI2337
    TOI-237 : TOI237
    TOI-2411 : TOI2411
    TOI-2427 : TOI2427
    TOI-257 : TOI257
    TOI-2669 : TOI2669
    TOI-269 : TOI269
    TOI-3362 : TOI3362
    TOI-421 : TOI421
    TOI-431 : TOI431
    TOI-4329 : TOI4329
    TOI-451 : TOI451
    TOI-481 : TOI481
    TOI-500 : TOI500
    TOI-530 : TOI530
    TOI-540 : TOI540
    TOI-544 : TOI544
    TOI-559 : TOI559
    TOI-561 : TOI561
    TOI-564 : TOI564
    TOI-620 : TOI620
    TOI-628 : TOI628
    TOI-640 : TOI640
    TOI-674 : TOI674
    TOI-677 : TOI677
    TOI-776 : TOI776
    TOI-813 : TOI813
    TOI-824 : TOI824
    TOI-833 : TOI833
    TOI-837 : TOI837
    TOI-892 : TOI892
    TOI-905 : TOI905
    TOI-954 : TOI954
    TrES-4 : TRES4
    TrES-5 : TRES5
    V1298 Tau : V1298TAU
    WASP-105 : WASP105
    WASP-106 : WASP106
    WASP-110 : WASP110
    WASP-113 : WASP113
    WASP-114 : WASP114
    WASP-117 : WASP117
    WASP-118 : WASP118
    WASP-119 : WASP119
    WASP-120 : WASP120
    WASP-123 : WASP123
    WASP-124 : WASP124
    WASP-126 : WASP126
    WASP-132 : WASP132
    WASP-133 : WASP133
    WASP-135 : WASP135
    WASP-136 : WASP136
    WASP-138 : WASP138
    WASP-139 : WASP139
    WASP-141 : WASP141
    WASP-142 : WASP142
    WASP-145 A : WASP145A
    WASP-147 : WASP147
    WASP-148 : WASP148
    WASP-151 : WASP151
    WASP-153 : WASP153
    WASP-156 : WASP156
    WASP-157 : WASP157
    WASP-158 : WASP158
    WASP-159 : WASP159
    WASP-160 B : WASP160B
    WASP-161 : WASP161
    WASP-163 : WASP163
    WASP-164 : WASP164
    WASP-165 : WASP165
    WASP-166 : WASP166
    WASP-167 : WASP167
    WASP-168 : WASP168
    WASP-169 : WASP169
    WASP-170 : WASP170
    WASP-172 : WASP172
    WASP-173 A : WASP173A
    WASP-174 : WASP174
    WASP-175 : WASP175
    WASP-176 : WASP176
    WASP-177 : WASP177
    WASP-178 : WASP178
    WASP-180 A : WASP180A
    WASP-181 : WASP181
    WASP-182 : WASP182
    WASP-183 : WASP183
    WASP-184 : WASP184
    WASP-185 : WASP185
    WASP-186 : WASP186
    WASP-187 : WASP187
    WASP-189 : WASP189
    WASP-190 : WASP190
    WASP-192 : WASP192
    WASP-20 : WASP20
    WASP-22 : WASP22
    WASP-23 : WASP23
    WASP-25 : WASP25
    WASP-41 : WASP41
    WASP-42 : WASP42
    WASP-44 : WASP44
    WASP-45 : WASP45
    WASP-47 : WASP47
    WASP-53 : WASP53
    WASP-54 : WASP54
    WASP-55 : WASP55
    WASP-56 : WASP56
    WASP-57 : WASP57
    WASP-58 : WASP58
    WASP-61 : WASP61
    WASP-66 : WASP66
    WASP-68 : WASP68
    WASP-70 A : WASP70A
    WASP-71 : WASP71
    WASP-73 : WASP73
    WASP-81 : WASP81
    WASP-82 : WASP82
    WASP-83 : WASP83
    WASP-84 : WASP84
    WASP-85 A : WASP85A
    WASP-88 : WASP88
    WASP-89 : WASP89
    WASP-90 : WASP90
    WASP-91 : WASP91
    WASP-92 : WASP92
    WASP-93 : WASP93
    WASP-98 : WASP98
    WASP-99 : WASP99
    Wolf 503 : WOLF503
    XO-6 : XO6
    XO-7 : XO7
    pi Men : PIMEN
    HATS-22 : HATS22
    K2-30 : K230
    Kepler-1308 : KEPLER1308
    Qatar-3 : QATAR3
    WASP-129 : WASP129
    WASP-144 : WASP144
    Kepler-51 : KEPLER51
    WD 1856+534 : WD1856534
    GJ 341 : GJ341
    GJ 4102 : GJ4102
    HD 80606 : HD80606
    GJ 4332 : GJ4332
    GJ 1008 : GJ1008
    LTT 5972 : LTT5972
    GJ 3090 : GJ3090
    GJ 806 : GJ806
    HD 109833 : HD109833
    HD 207496 : HD207496
    HD 260655 : HD260655
    HD 93963 A : HD93963A
    HD 95338 : HD95338
    HIP 29442 : HIP29442
    HIP 94235 : HIP94235
    HIP 97166 : HIP97166
    K2-105 : K2105
    K2-370 : K237
    K2-415 : K2415
    K2-60 : K260
    Kepler-1656 : Kepler1656
    Kepler-63 : Kepler63
    Kepler-96 : Kepler96
    NGTS-14 A : NGTS14A
    TOI-1136 : TOI1136
    TOI-1272 : TOI1272
    TOI-1278 : TOI1278
    TOI-1288 : TOI1288
    TOI-132 : TOI132
    TOI-139 : TOI139
    TOI-1422 : TOI1422
    TOI-1468 : TOI1468
    TOI-1470 : TOI1470
    TOI-1634 : TOI1634
    TOI-1694 : TOI1694
    TOI-1695 : TOI1695
    TOI-1710 : TOI1710
    TOI-1801 : TOI1801
    TOI-181 : TOI181
    TOI-1853 : TOI1853
    TOI-1859 : TOI1859
    TOI-199 : TOI199
    TOI-2000 : TOI2000
    TOI-2010 : TOI2010
    TOI-2018 : TOI2018
    TOI-2134 : TOI2134
    TOI-2136 : TOI2136
    TOI-220 : TOI220
    TOI-2364 : TOI2364
    TOI-2443 : TOI2443
    TOI-2445 : TOI2445
    TOI-2459 : TOI2459
    TOI-2498 : TOI2498
    TOI-251 : TOI251
    TOI-277 : TOI277
    TOI-3082 : TOI3082
    TOI-332 : TOI332
    TOI-3629 : TOI3629
    TOI-3785 : TOI3785
    TOI-4010 : TOI4010
    TOI-444 : TOI444
    TOI-4479 : TOI4479
    TOI-4600 : TOI4600
    TOI-4641 : TOI4641
    TOI-470 : TOI470
    TOI-5126 : TOI5126
    TOI-532 : TOI532
    TOI-5344 : TOI5344
    TOI-5398 : TOI5398
    TOI-5678 : TOI5678
    TOI-5704 : TOI5704
    TOI-5803 : TOI5803
    TOI-672 : TOI672
    TOI-712 : TOI712
    TOI-908 : TOI908
    TOI-913 : TOI913
    TOI-942 : TOI942
    TOI-969 : TOI969
    Wolf 327 : Wolf327
    CoRoT-3 : COROT3
    CoRoT-36 : COROT36
    Gaia-1 : GAIA1
    Gaia-2 : GAIA2
    HATS-10 : HATS10
    HATS-12 : HATS12
    HATS-45 : HATS45
    HATS-55 : HATS55
    HATS-61 : HATS61
    HD 118203 : HD118203
    HD 15906 : HD15906
    HD 21749 : HD21749
    HD 235088 : HD235088
    HD 28109 : HD28109
    HIP 113103 : HIP113103
    HIP 116454 : HIP116454
    HIP 9618 : HIP9618
    K2-233 : K2233
    K2-240 : K2240
    K2-285 : K2285
    K2-321 : K2321
    K2-344 : K2344
    K2-348 : K2348
    K2-417 : K2417
    K2-99 : K299
    Kepler-1515 : KEPLER1515
    Kepler-1517 : KEPLER1517
    Kepler-1658 : KEPLER1658
    Kepler-411 : KEPLER411
    Kepler-91 : KEPLER91
    KOI-12 : KOI12
    NGTS-24 : NGTS24
    NGTS-9 : NGTS9
    PH2 : PH2
    TIC 237913194 : TIC237913194
    TOI-1107 : TOI1107
    TOI-1181 : TOI1181
    TOI-1194 : TOI1194
    TOI-1246 : TOI1246
    TOI-1338 : TOI1338
    TOI-1408 : TOI1408
    TOI-1416 : TOI1416
    TOI-1420 : TOI1420
    TOI-1452 : TOI1452
    TOI-1516 : TOI1516
    TOI-1680 : TOI1680
    TOI-1811 : TOI1811
    TOI-1820 : TOI1820
    TOI-1937 A : TOI1937A
    TOI-198 : TOI198
    TOI-2025 : TOI2025
    TOI-2046 : TOI2046
    TOI-2048 : TOI2048
    TOI-206 : TOI206
    TOI-2081 : TOI2081
    TOI-2145 : TOI2145
    TOI-2152 A : TOI2152A
    TOI-2154 : TOI2154
    TOI-2158 : TOI2158
    TOI-2193 A : TOI2193A
    TOI-2194 : TOI2194
    TOI-2202 : TOI2202
    TOI-2207 : TOI2207
    TOI-2236 : TOI2236
    TOI-2338 : TOI2338
    TOI-2421 : TOI2421
    TOI-2497 : TOI2497
    TOI-2524 : TOI2524
    TOI-2567 : TOI2567
    TOI-2570 : TOI2570
    TOI-2583 A : TOI2583A
    TOI-2587 A : TOI2587A
    TOI-262 : TOI262
    TOI-2641 : TOI2641
    TOI-2796 : TOI2796
    TOI-2803 A : TOI2803A
    TOI-2818 : TOI2818
    TOI-2842 : TOI2842
    TOI-2977 : TOI2977
    TOI-3023 : TOI3023
    TOI-3235 : TOI3235
    TOI-3331 A : TOI3331A
    TOI-3364 : TOI3364
    TOI-3540 A : TOI3540A
    TOI-3688 A : TOI3688A
    TOI-3693 : TOI3693
    TOI-3714 : TOI3714
    TOI-3807 : TOI3807
    TOI-3819 : TOI3819
    TOI-3884 : TOI3884
    TOI-3912 : TOI3912
    TOI-3976 A : TOI3976A
    TOI-4087 : TOI4087
    TOI-411 : TOI411
    TOI-4137 : TOI4137
    TOI-4145 A : TOI4145A
    TOI-4308 : TOI4308
    TOI-4342 : TOI4342
    TOI-4377 : TOI4377
    TOI-4406 : TOI4406
    TOI-4463 A : TOI4463A
    TOI-4551 : TOI4551
    TOI-4603 : TOI4603
    TOI-4791 : TOI4791
    TOI-615 : TOI615
    TOI-622 : TOI622
    TOI-778 : TOI778
    TOI-858 : TOI858
    WASP-130 : WASP130
    WASP-150 : WASP150
    WASP-162 : WASP162
    WASP-171 : WASP171
    WASP-193 : WASP193
    WASP-60 : WASP60
    HD 110067 : HD110067
    TOI-3261 : TOI3261
    CoRoT-10 : CoRoT10
    CoRoT-12 : CoRoT12
    CoRoT-13 : CoRoT13
    CoRoT-14 : CoRoT14
    CoRoT-16 : CoRoT16
    CoRoT-18 : CoRoT18
    CoRoT-22 : CoRoT22
    CoRoT-24 : CoRoT24
    CoRoT-25 : CoRoT25
    CoRoT-26 : CoRoT26
    CoRoT-28 : CoRoT28
    CoRoT-29 : CoRoT29
    CoRoT-31 : CoRoT31
    CoRoT-4 : CoRoT4
    CoRoT-6 : CoRoT6
    CoRoT-8 : CoRoT8
    CoRoT-9 : CoRoT9
    CoRoTID 223977153 : CoRoTID223977153
    EPIC 212297394 : EPIC212297394
    EPIC 220674823 : EPIC220674823
    EPIC 229004835 : EPIC229004835
    EPIC 249893012 : EPIC249893012
    HATS-14 : HATS14
    HATS-15 : HATS15
    HATS-16 : HATS16
    HATS-17 : HATS17
    HATS-32 : HATS32
    HATS-36 : HATS36
    HATS-44 : HATS44
    HATS-49 : HATS49
    HATS-54 : HATS54
    HATS-59 : HATS59
    HATS-63 : HATS63
    HATS-69 : HATS69
    HATS-71 : HATS71
    HATS-74 A : HATS74A
    HATS-75 : HATS75
    HATS-76 : HATS76
    HATS-77 : HATS77
    HATS-8 : HATS8
    HD 137496 : HD137496
    HD 18599 : HD18599
    HD 20329 : HD20329
    HD 80653 : HD80653
    HD 88986 : HD88986
    K2-10 : K210
    K2-100 : K2100
    K2-108 : K2108
    K2-110 : K2110
    K2-111 : K2111
    K2-113 : K2113
    K2-114 : K2114
    K2-115 : K2115
    K2-117 : K2117
    K2-12 : K212
    K2-122 : K2122
    K2-126 : K2126
    K2-127 : K2127
    K2-128 : K2128
    K2-130 : K2130
    K2-131 : K2131
    K2-133 : K2133
    K2-14 : K214
    K2-146 : K2146
    K2-148 : K2148
    K2-151 : K2151
    K2-152 : K2152
    K2-153 : K2153
    K2-154 : K2154
    K2-158 : K2158
    K2-159 : K2159
    K2-160 : K2160
    K2-161 : K2161
    K2-162 : K2162
    K2-163 : K2163
    K2-164 : K2164
    K2-165 : K2165
    K2-168 : K2168
    K2-17 : K217
    K2-170 : K2170
    K2-171 : K2171
    K2-172 : K2172
    K2-173 : K2173
    K2-175 : K2175
    K2-176 : K2176
    K2-177 : K2177
    K2-178 : K2178
    K2-179 : K2179
    K2-180 : K2180
    K2-181 : K2181
    K2-182 : K2182
    K2-184 : K2184
    K2-186 : K2186
    K2-188 : K2188
    K2-189 : K2189
    K2-190 : K2190
    K2-193 : K2193
    K2-194 : K2194
    K2-195 : K2195
    K2-196 : K2196
    K2-197 : K2197
    K2-199 : K2199
    K2-200 : K2200
    K2-201 : K2201
    K2-203 : K2203
    K2-204 : K2204
    K2-206 : K2206
    K2-208 : K2208
    K2-216 : K2216
    K2-216 : K2216
    K2-217 : K2217
    K2-219 : K2219
    K2-220 : K2220
    K2-221 : K2221
    K2-224 : K2224
    K2-225 : K2225
    K2-226 : K2226
    K2-227 : K2227
    K2-228 : K2228
    K2-229 : K2229
    K2-230 : K2230
    K2-241 : K2241
    K2-242 : K2242
    K2-245 : K2245
    K2-246 : K2246
    K2-247 : K2247
    K2-248 : K2248
    K2-250 : K2250
    K2-253 : K2253
    K2-254 : K2254
    K2-255 : K2255
    K2-265 : K2265
    K2-268 : K2268
    K2-27 : K227
    K2-270 : K2270
    K2-271 : K2271
    K2-272 : K2272
    K2-273 : K2273
    K2-274 : K2274
    K2-275 : K2275
    K2-276 : K2276
    K2-277 : K2277
    K2-281 : K2281
    K2-283 : K2283
    K2-286 : K2286
    K2-291 : K2291
    K2-292 : K2292
    K2-294 : K2294
    K2-308 : K2308
    K2-331 : K2331
    K2-341 : K2341
    K2-342 : K2342
    K2-346 : K2346
    K2-35 : K235
    K2-352 : K2352
    K2-354 : K2354
    K2-357 : K2357
    K2-365 : K2365
    K2-366 : K2366
    K2-367 : K2367
    K2-368 : K2368
    K2-369 : K2369
    K2-37 : K237
    K2-371 : K2371
    K2-372 : K2372
    K2-374 : K2374
    K2-376 : K2376
    K2-378 : K2378
    K2-379 : K2379
    K2-38 : K238
    K2-380 : K2380
    K2-381 : K2381
    K2-383 : K2383
    K2-387 : K2387
    K2-388 : K2388
    K2-389 : K2389
    K2-390 : K2390
    K2-391 : K2391
    K2-393 : K2393
    K2-394 : K2394
    K2-395 : K2395
    K2-396 : K2396
    K2-397 : K2397
    K2-398 : K2398
    K2-401 : K2401
    K2-404 : K2404
    K2-407 : K2407
    K2-408 : K2408
    K2-409 : K2409
    K2-416 : K2416
    K2-43 : K243
    K2-44 : K244
    K2-45 : K245
    K2-5 : K25
    K2-62 : K262
    K2-65 : K265
    K2-68 : K268
    K2-70 : K270
    K2-73 : K273
    K2-75 : K275
    K2-77 : K277
    K2-8 : K28
    K2-80 : K280
    K2-81 : K281
    K2-84 : K284
    K2-86 : K286
    K2-98 : K298
    KIC 9663113 : KIC9663113
    KOI-1257 : KOI1257
    KOI-217 : KOI217
    KOI-351 : KOI351
    KOI-3680 : KOI3680
    KOI-7368 : KOI7368
    KOI-984 : KOI984
    Kepler-1002 : Kepler1002
    Kepler-1003 : Kepler1003
    Kepler-1004 : Kepler1004
    Kepler-101 : Kepler101
    Kepler-1011 : Kepler1011
    Kepler-1012 : Kepler1012
    Kepler-1016 : Kepler1016
    Kepler-1019 : Kepler1019
    Kepler-1028 : Kepler1028
    Kepler-103 : Kepler103
    Kepler-1035 : Kepler1035
    Kepler-1037 : Kepler1037
    Kepler-1039 : Kepler1039
    Kepler-1055 : Kepler1055
    Kepler-1061 : Kepler1061
    Kepler-1065 : Kepler1065
    Kepler-1066 : Kepler1066
    Kepler-107 : Kepler107
    Kepler-1072 : Kepler1072
    Kepler-1073 : Kepler1073
    Kepler-1075 : Kepler1075
    Kepler-1082 : Kepler1082
    Kepler-109 : Kepler109
    Kepler-1091 : Kepler1091
    Kepler-1099 : Kepler1099
    Kepler-1106 : Kepler1106
    Kepler-1107 : Kepler1107
    Kepler-1139 : Kepler1139
    Kepler-1148 : Kepler1148
    Kepler-1167 : Kepler1167
    Kepler-117 : Kepler117
    Kepler-1179 : Kepler1179
    Kepler-1181 : Kepler1181
    Kepler-1189 : Kepler1189
    Kepler-1193 : Kepler1193
    Kepler-1205 : Kepler1205
    Kepler-1206 : Kepler1206
    Kepler-1217 : Kepler1217
    Kepler-1228 : Kepler1228
    Kepler-1256 : Kepler1256
    Kepler-1259 : Kepler1259
    Kepler-1270 : Kepler1270
    Kepler-1271 : Kepler1271
    Kepler-1310 : Kepler1310
    Kepler-1311 : Kepler1311
    Kepler-1312 : Kepler1312
    Kepler-1313 : Kepler1313
    Kepler-1315 : Kepler1315
    Kepler-1317 : Kepler1317
    Kepler-1320 : Kepler1320
    Kepler-1323 : Kepler1323
    Kepler-1344 : Kepler1344
    Kepler-1346 : Kepler1346
    Kepler-1356 : Kepler1356
    Kepler-148 : Kepler148
    Kepler-15 : Kepler15
    Kepler-1514 : Kepler1514
    Kepler-1519 : Kepler1519
    Kepler-1520 : Kepler1520
    Kepler-1522 : Kepler1522
    Kepler-1523 : Kepler1523
    Kepler-1524 : Kepler1524
    Kepler-1528 : Kepler1528
    Kepler-1530 : Kepler1530
    Kepler-154 : Kepler154
    Kepler-157 : Kepler157
    Kepler-1604 : Kepler1604
    Kepler-1613 : Kepler1613
    Kepler-1624 : Kepler1624
    Kepler-1627 : Kepler1627
    Kepler-1642 : Kepler1642
    Kepler-1654 : Kepler1654
    Kepler-1655 : Kepler1655
    Kepler-166 : Kepler166
    Kepler-1704 : Kepler1704
    Kepler-198 : Kepler198
    Kepler-289 : Kepler289
    Kepler-30 : Kepler30
    Kepler-323 : Kepler323
    Kepler-4 : Kepler4
    Kepler-41 : Kepler41
    Kepler-412 : Kepler412
    Kepler-418 : Kepler418
    Kepler-421 : Kepler421
    Kepler-422 : Kepler422
    Kepler-423 : Kepler423
    Kepler-424 : Kepler424
    Kepler-425 : Kepler425
    Kepler-426 : Kepler426
    Kepler-427 : Kepler427
    Kepler-428 : Kepler428
    Kepler-43 : Kepler43
    Kepler-44 : Kepler44
    Kepler-46 : Kepler46
    Kepler-460 : Kepler460
    Kepler-463 : Kepler463
    Kepler-464 : Kepler464
    Kepler-466 : Kepler466
    Kepler-472 : Kepler472
    Kepler-473 : Kepler473
    Kepler-474 : Kepler474
    Kepler-475 : Kepler475
    Kepler-476 : Kepler476
    Kepler-485 : Kepler485
    Kepler-487 : Kepler487
    Kepler-490 : Kepler490
    Kepler-491 : Kepler491
    Kepler-495 : Kepler495
    Kepler-497 : Kepler497
    Kepler-498 : Kepler498
    Kepler-501 : Kepler501
    Kepler-502 : Kepler502
    Kepler-507 : Kepler507
    Kepler-511 : Kepler511
    Kepler-514 : Kepler514
    Kepler-518 : Kepler518
    Kepler-523 : Kepler523
    Kepler-529 : Kepler529
    Kepler-530 : Kepler530
    Kepler-531 : Kepler531
    Kepler-532 : Kepler532
    Kepler-533 : Kepler533
    Kepler-535 : Kepler535
    Kepler-536 : Kepler536
    Kepler-537 : Kepler537
    Kepler-539 : Kepler539
    Kepler-541 : Kepler541
    Kepler-543 : Kepler543
    Kepler-546 : Kepler546
    Kepler-547 : Kepler547
    Kepler-548 : Kepler548
    Kepler-550 : Kepler550
    Kepler-552 : Kepler552
    Kepler-553 : Kepler553
    Kepler-554 : Kepler554
    Kepler-557 : Kepler557
    Kepler-56 : Kepler56
    Kepler-561 : Kepler561
    Kepler-562 : Kepler562
    Kepler-564 : Kepler564
    Kepler-565 : Kepler565
    Kepler-568 : Kepler568
    Kepler-578 : Kepler578
    Kepler-585 : Kepler585
    Kepler-586 : Kepler586
    Kepler-590 : Kepler590
    Kepler-592 : Kepler592
    Kepler-605 : Kepler605
    Kepler-608 : Kepler608
    Kepler-609 : Kepler609
    Kepler-611 : Kepler611
    Kepler-612 : Kepler612
    Kepler-618 : Kepler618
    Kepler-619 : Kepler619
    Kepler-620 : Kepler620
    Kepler-621 : Kepler621
    Kepler-629 : Kepler629
    Kepler-634 : Kepler634
    Kepler-636 : Kepler636
    Kepler-645 : Kepler645
    Kepler-648 : Kepler648
    Kepler-650 : Kepler650
    Kepler-656 : Kepler656
    Kepler-664 : Kepler664
    Kepler-666 : Kepler666
    Kepler-667 : Kepler667
    Kepler-669 : Kepler669
    Kepler-670 : Kepler670
    Kepler-673 : Kepler673
    Kepler-674 : Kepler674
    Kepler-677 : Kepler677
    Kepler-678 : Kepler678
    Kepler-680 : Kepler680
    Kepler-682 : Kepler682
    Kepler-683 : Kepler683
    Kepler-684 : Kepler684
    Kepler-685 : Kepler685
    Kepler-686 : Kepler686
    Kepler-688 : Kepler688
    Kepler-690 : Kepler690
    Kepler-693 : Kepler693
    Kepler-694 : Kepler694
    Kepler-695 : Kepler695
    Kepler-696 : Kepler696
    Kepler-697 : Kepler697
    Kepler-700 : Kepler700
    Kepler-702 : Kepler702
    Kepler-703 : Kepler703
    Kepler-707 : Kepler707
    Kepler-708 : Kepler708
    Kepler-714 : Kepler714
    Kepler-715 : Kepler715
    Kepler-718 : Kepler718
    Kepler-719 : Kepler719
    Kepler-720 : Kepler720
    Kepler-722 : Kepler722
    Kepler-723 : Kepler723
    Kepler-724 : Kepler724
    Kepler-725 : Kepler725
    Kepler-728 : Kepler728
    Kepler-729 : Kepler729
    Kepler-731 : Kepler731
    Kepler-734 : Kepler734
    Kepler-736 : Kepler736
    Kepler-74 : Kepler74
    Kepler-740 : Kepler740
    Kepler-743 : Kepler743
    Kepler-75 : Kepler75
    Kepler-750 : Kepler750
    Kepler-755 : Kepler755
    Kepler-762 : Kepler762
    Kepler-767 : Kepler767
    Kepler-77 : Kepler77
    Kepler-773 : Kepler773
    Kepler-775 : Kepler775
    Kepler-78 : Kepler78
    Kepler-782 : Kepler782
    Kepler-785 : Kepler785
    Kepler-796 : Kepler796
    Kepler-799 : Kepler799
    Kepler-80 : Kepler80
    Kepler-805 : Kepler805
    Kepler-808 : Kepler808
    Kepler-814 : Kepler814
    Kepler-815 : Kepler815
    Kepler-816 : Kepler816
    Kepler-817 : Kepler817
    Kepler-818 : Kepler818
    Kepler-822 : Kepler822
    Kepler-825 : Kepler825
    Kepler-827 : Kepler827
    Kepler-828 : Kepler828
    Kepler-831 : Kepler831
    Kepler-842 : Kepler842
    Kepler-843 : Kepler843
    Kepler-845 : Kepler845
    Kepler-849 : Kepler849
    Kepler-855 : Kepler855
    Kepler-856 : Kepler856
    Kepler-857 : Kepler857
    Kepler-858 : Kepler858
    Kepler-860 : Kepler860
    Kepler-867 : Kepler867
    Kepler-872 : Kepler872
    Kepler-890 : Kepler890
    Kepler-891 : Kepler891
    Kepler-893 : Kepler893
    Kepler-904 : Kepler904
    Kepler-905 : Kepler905
    Kepler-908 : Kepler908
    Kepler-912 : Kepler912
    Kepler-914 : Kepler914
    Kepler-915 : Kepler915
    Kepler-922 : Kepler922
    Kepler-932 : Kepler932
    Kepler-943 : Kepler943
    Kepler-950 : Kepler950
    Kepler-951 : Kepler951
    Kepler-952 : Kepler952
    Kepler-953 : Kepler953
    Kepler-954 : Kepler954
    Kepler-956 : Kepler956
    Kepler-957 : Kepler957
    Kepler-960 : Kepler960
    Kepler-972 : Kepler972
    Kepler-975 : Kepler975
    Kepler-990 : Kepler990
    Kepler-994 : Kepler994
    Kepler-996 : Kepler996
    Kepler-997 : Kepler997
    LHS 1815 : LHS1815
    NGTS-1 : NGTS1
    NGTS-15 : NGTS15
    NGTS-16 : NGTS16
    NGTS-17 : NGTS17
    NGTS-18 : NGTS18
    NGTS-20 : NGTS20
    NGTS-21 : NGTS21
    NGTS-23 : NGTS23
    NGTS-25 : NGTS25
    NGTS-3 A : NGTS3A
    NGTS-4 : NGTS4
    OGLE-TR-182 : OGLETR182
    POTS-1 : POTS1
    TIC 172900988 : TIC172900988
    TIC 279401253 : TIC279401253
    TOI-1052 : TOI1052
    TOI-1062 : TOI1062
    TOI-1221 : TOI1221
    TOI-1235 : TOI1235
    TOI-1444 : TOI1444
    TOI-1696 : TOI1696
    TOI-1736 : TOI1736
    TOI-2084 : TOI2084
    TOI-2095 : TOI2095
    TOI-2096 : TOI2096
    TOI-2141 : TOI2141
    TOI-2180 : TOI2180
    TOI-2184 : TOI2184
    TOI-2196 : TOI2196
    TOI-2257 : TOI2257
    TOI-2285 : TOI2285
    TOI-244 : TOI244
    TOI-2525 : TOI2525
    TOI-2589 : TOI2589
    TOI-3757 : TOI3757
    TOI-3984 A : TOI3984A
    TOI-4127 : TOI4127
    TOI-4184 : TOI4184
    TOI-4201 : TOI4201
    TOI-4562 : TOI4562
    TOI-4582 : TOI4582
    TOI-4860 : TOI4860
    TOI-5205 : TOI5205
    TOI-5293 A : TOI5293A
    TOI-5542 : TOI5542
    TOI-733 : TOI733
    TOI-763 : TOI763
    TOI-784 : TOI784
    WASP-59 : WASP59
    WTS-2 : WTS2
    Wendelstein-1 : Wendelstein1
    Wendelstein-2 : Wendelstein2
    TOI-6086 : TOI6086
    GJ 12 : GJ12
    TOI-1135 : TOI1135
    TOI-4336 A : TOI336A
    SPECULOOS-3 : SPECULOOS3
    TOI-904 : TOI904
    TOI-4559 : TOI4559
    TOI-771 : TOI771
    TOI-6008 : TOI6008
    TOI-5799 : TOI5799
    TOI-2120 : TOI2120
    TOI-782 : TOI782
    TOI-260 : TOI260
    '''
    return


# -------------------- -----------------------------------------------
