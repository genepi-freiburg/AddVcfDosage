#!/bin/bash
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java \
	-cp /data/gwas/scripts/AddVcfDosage/lib/htsjdk-2.5.0.jar:/data/gwas/scripts/AddVcfDosage/lib/AddVcfDosage.jar \
	 de.wuttke.vcf.AddVcfDosage "$@"
