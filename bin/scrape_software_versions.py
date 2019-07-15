#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'UCL-BLIC/rnaseq_variant_calling': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'STAR': ['v_star.txt', r"(\S+)"],
    'Picard MarkDuplicates': ['v_markduplicates.txt', r"([\d\.]+)-SNAPSHOT"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['UCL-BLIC/rnaseq_variant_calling'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = False
results['Picard MarkDuplicates'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = match.group(1)

# Strip STAR or HiSAT2
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'RNAseq Variant Calling Software Versions'
section_href: 'https://github.com/UCL-BLIC/rnaseq_variant_calling'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")
