"""This is a temporary module to synchronise the parameters/options lists with IMPRESS.
Given the size of the files (CategoryRemapping.tsv, MergeParameterList.txt), it's impossible
to modify them by hand, so this module is used instead. In the future, parameter mapping
and representation will be improved and refactored, so this module will no longer be necessary."""

import csv

print("Load CategoryMap.conf")
category_map: dict = dict()
for line in open(
    "Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/CategoryMap.conf"
):
    line = line.rstrip()
    if ":" in line:
        a, b = line.split(":")
        category_map[a] = b

print("Load MergeParameterList.txt")
merge_parameter_list = set(
    open(
        "Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/MergeParameterList.txt"
    )
    .read()
    .splitlines()[1:]
)

print("Load CategoryRemapping.tsv")
# Center, parameter, value: assignment.
c_p_v_a: dict = {}
with open(
    "Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/annotation/CategoryRemapping.tsv",
    mode="r",
    newline="",
) as file:
    file.readline()
    for row in file:
        centre_parameter, value, assignment = row.rstrip().split("\t")
        if not "_EYE_" in centre_parameter:
            continue
        parts = centre_parameter.split("_")
        centre = parts[0]
        parameter = "_".join(parts[1:])
        c_p_v_a[f"{centre}/{parameter}/{value}"] = int(assignment)
# Parameter, value: assignment.
p_v_a: dict = {}
for c_p_v, a in c_p_v_a.items():
    c, p, v = c_p_v.split("/")
    key = f"{p}_{v}"
    if key in p_v_a:
        assert p_v_a[key] == a
    else:
        p_v_a[key] = a
# Add manual new assignments.
p_v_a["EYE_092_002_normal"] = 0
p_v_a["EYE_092_002_both eyes abnormal"] = 1
p_v_a["EYE_092_002_right eye abnormal"] = 1
p_v_a["EYE_092_002_left eye abnormal"] = 1

print("Process the new file")
merge_parameter_list_out = open(
    "Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/MergeParameterList.txt",
    "a",
)
category_remapping_out = open(
    "Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/annotation/CategoryRemapping.tsv",
    "a",
)
with open("LA_EYE_parameters.csv", mode="r", newline="") as file:
    reader = csv.DictReader(file)
    for row in reader:
        # Extract centre and the actual parameter.
        centre_parameter = row["parameter_stable_id"]
        parts = centre_parameter.split("_")
        centre = parts[0]
        parameter = "_".join(parts[1:])
        # Amend MergeParameterList, if needed.
        if centre_parameter not in merge_parameter_list:
            merge_parameter_list_out.write(centre_parameter + "\n")
        # Amend CategoryRemapping.tsv, if needed.
        for value in sorted(row["options"].split(";")):
            renamed = category_map.get(value, value)
            if renamed == "NA":
                continue
            key = f"{centre}/{parameter}/{renamed}"
            if not key in c_p_v_a:
                a = p_v_a[f"{parameter}_{renamed}"]
                category_remapping_out.write(f"{centre_parameter}\t{renamed}\t{a}\n")
                c_p_v_a[key] = a

merge_parameter_list_out.close()
category_remapping_out.close()
