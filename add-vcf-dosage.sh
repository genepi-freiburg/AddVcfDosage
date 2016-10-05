#!/bin/bash
JAVA=/data/ngs/bin/jdk1.8.0_101/bin/java
SCRIPT_LIB=/data/gwas/scripts/AddVcfDosage/lib
$JAVA \
	-cp $SCRIPT_LIB/htsjdk-2.5.0.jar:$SCRIPT_LIB/AddVcfDosage.jar \
	 de.wuttke.vcf.AddVcfDosage "$@"
