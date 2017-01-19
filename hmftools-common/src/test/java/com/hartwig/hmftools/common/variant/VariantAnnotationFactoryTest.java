package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VariantAnnotationFactoryTest {

    @Test
    public void canLoadFromEmptyVCFInfoField() {
        final List<VariantAnnotation> annotations = VariantAnnotationFactory.fromVCFInfoField(Strings.EMPTY);
        assertEquals(0, annotations.size());
    }

    @Test
    public void canDealWithMultipleConsequences() {
        final String info = "ANN=allele|missense_variant&splice_donor_variant";
        final List<VariantAnnotation> annotations = VariantAnnotationFactory.fromVCFInfoField(info);
        assertEquals(1, annotations.size());
        assertEquals(VariantConsequence.MISSENSE_VARIANT, annotations.get(0).consequences().get(0));
        assertEquals(VariantConsequence.SPLICE_DONOR_VARIANT, annotations.get(0).consequences().get(1));
    }

    @Test
    public void canLoadFromTrivialVCFInfoField() {
        final String info = "ANN=allele|consequence|severity|gene|geneID|featureType|featureID|"
                + "transcriptBioType|rank|hgvsCoding|hgvsProtein|cDNAPosAndLength|cdsPosAndLength|"
                + "aaPosAndLength|distance|addition";
        final List<VariantAnnotation> annotations = VariantAnnotationFactory.fromVCFInfoField(info);
        assertEquals(1, annotations.size());
        final VariantAnnotation annotation = annotations.get(0);
        assertEquals("allele", annotation.allele());
        assertEquals(VariantConsequence.OTHER, annotation.consequences().get(0));
        assertEquals("severity", annotation.severity());
        assertEquals("gene", annotation.gene());
        assertEquals("geneID", annotation.geneID());
        assertEquals("featureType", annotation.featureType());
        assertEquals("featureID", annotation.featureID());
        assertEquals("transcriptBioType", annotation.transcriptBioType());
        assertEquals("rank", annotation.rank());
        assertEquals("hgvsCoding", annotation.hgvsCoding());
        assertEquals("hgvsProtein", annotation.hgvsProtein());
        assertEquals("cDNAPosAndLength", annotation.cDNAPosAndLength());
        assertEquals("cdsPosAndLength", annotation.cdsPosAndLength());
        assertEquals("aaPosAndLength", annotation.aaPosAndLength());
        assertEquals("distance", annotation.distance());
        assertEquals("addition", annotation.addition());
    }

    @Test
    public void canLoadFromRealVCFInfoField() {
        final String info = "AB=0.119266;ABP=140.251;AC=2;AF=0.286;AN=7;"
                + "ANN=C|missense_variant|MODERATE|NOTCH1|ENSG00000148400|transcript|ENST00000277541|protein_coding|"
                + "18/34|c.2915C>G|p.Pro972Arg|2991/9371|2915/7668|972/2555||,C|sequence_feature|LOW|NOTCH1|"
                + "ENSG00000148400|topological_domain:Extracellular|ENST00000277541|protein_coding||c.2915C>G||||||,"
                + "C|sequence_feature|LOW|NOTCH1|ENSG00000148400|domain:EGF-like_25|ENST00000277541|protein_coding||"
                + "c.2915C>G||||||,C|sequence_feature|LOW|NOTCH1|ENSG00000148400|disulfide_bond|ENST00000277541|"
                + "protein_coding||c.2915C>G||||||;AO=13;CIGAR=1X;DP=124;DPB=124;DPRA=7.26667;EPP=3.17734;"
                + "EPPR=3.96888;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=2;NUMALT=1;ODDS=11.378;PAIRED=1;PAIREDR=1;"
                + "PAO=0;PQA=0;PQR=0;PRO=0;QA=482;QR=4272;RO=111;RPL=9;RPP=7.18621;RPPR=3.18637;RPR=4;RUN=1;SAF=6;"
                + "SAP=3.17734;SAR=7;SOMATIC;SRF=55;SRP=3.02986;SRR=56;TYPE=snp;VT=SNP;set=freebayes-mutect;"
                + "technology.ILLUMINA=1;CSA=2,2;CSP=2";
        final List<VariantAnnotation> annotations = VariantAnnotationFactory.fromVCFInfoField(info);
        assertEquals(4, annotations.size());
        assertEquals(VariantConsequence.MISSENSE_VARIANT, annotations.get(0).consequences().get(0));
    }
}