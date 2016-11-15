#!/bin/bash
JAVA=/usr/bin/java
SCRIPT_LIB=./lib
$JAVA \
	-cp $SCRIPT_LIB/htsjdk-2.5.0.jar:$SCRIPT_LIB/AddVcfDosage.jar \
	 de.wuttke.vcf.AddVcfDosage "$@"
