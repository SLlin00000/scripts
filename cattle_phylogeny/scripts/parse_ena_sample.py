#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as ET

xml_file = sys.argv[1]

def clean(x):
    if x is None:
        return ""
    return " ".join(x.strip().split())

try:
    root = ET.parse(xml_file).getroot()
except Exception:
    print("\t".join([""] * 9))
    sys.exit(0)

sample = root.find(".//SAMPLE")
if sample is None:
    sample = root

scientific_name = clean(sample.findtext(".//SCIENTIFIC_NAME"))
description = clean(sample.findtext(".//DESCRIPTION"))

attrs = {}

for attr in sample.findall(".//SAMPLE_ATTRIBUTE"):
    tag = clean(attr.findtext("TAG"))
    value = clean(attr.findtext("VALUE"))
    if tag:
        attrs[tag.lower()] = value

def get_attr(*names):
    for name in names:
        value = attrs.get(name.lower(), "")
        if value:
            return value
    return ""

isolate = get_attr("isolate")
breed = get_attr("breed", "breed name", "breed_name")
sex = get_attr("sex", "gender")
country = get_attr(
    "country",
    "geographic location (country and/or sea)",
    "geo_loc_name"
)
sample_name = get_attr(
    "sample name",
    "sample_name",
    "individual",
    "animal id",
    "animal_id"
)
dev_stage = get_attr(
    "dev_stage",
    "developmental stage",
    "development stage"
)

organism_attr = get_attr("organism")

print("\t".join([
    scientific_name,
    organism_attr,
    isolate,
    breed,
    sex,
    country,
    sample_name,
    dev_stage,
    description
]))

