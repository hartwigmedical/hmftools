package com.hartwig.hmftools.common.peach;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PeachGenotypeFileTest
{
    private static final String PEACH_PYTHON_GENOTYPES_FILE = Resources.getResource("peach/python/sample_genotype.tsv").getPath();
    private static final String PEACH_JAVA_GENOTYPES_FILE = Resources.getResource("peach/java/sample_genotype.tsv").getPath();

    @Test
    public void loadPeachPythonGenotypeFile() throws IOException
    {
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(PEACH_PYTHON_GENOTYPES_FILE);

        assertEquals(3, peachGenotypes.size());

        assertEquals("GENE", peachGenotypes.get(0).gene());
        assertEquals("*1_HET", peachGenotypes.get(0).haplotype());
        assertEquals("*1", peachGenotypes.get(0).allele());
        assertEquals(1, peachGenotypes.get(0).alleleCount());
        assertEquals("Reduced Function", peachGenotypes.get(0).function());
        assertEquals("Capecitabine;otherDrug", peachGenotypes.get(0).linkedDrugs());
        assertEquals("link;otherLink", peachGenotypes.get(0).urlPrescriptionInfo());
        assertEquals("panel", peachGenotypes.get(0).panelVersion());
        assertEquals("pilot", peachGenotypes.get(0).repoVersion());

        assertEquals("GENE2", peachGenotypes.get(1).gene());
        assertEquals("*3_HOM", peachGenotypes.get(1).haplotype());
        assertEquals("*3", peachGenotypes.get(1).allele());
        assertEquals(2, peachGenotypes.get(1).alleleCount());
        assertEquals("Normal Function", peachGenotypes.get(1).function());
        assertEquals("Capecitabine", peachGenotypes.get(1).linkedDrugs());
        assertEquals("link", peachGenotypes.get(1).urlPrescriptionInfo());
        assertEquals("panel", peachGenotypes.get(1).panelVersion());
        assertEquals("pilot", peachGenotypes.get(1).repoVersion());

        assertEquals("GENE3", peachGenotypes.get(2).gene());
        assertEquals("Unresolved Haplotype", peachGenotypes.get(2).haplotype());
        assertEquals("Unresolved Haplotype", peachGenotypes.get(2).allele());
        assertEquals(2, peachGenotypes.get(2).alleleCount());
        assertEquals("Unknown Function", peachGenotypes.get(2).function());
        assertEquals("Capecitabine", peachGenotypes.get(2).linkedDrugs());
        assertEquals("link", peachGenotypes.get(2).urlPrescriptionInfo());
        assertEquals("panel", peachGenotypes.get(2).panelVersion());
        assertEquals("pilot", peachGenotypes.get(2).repoVersion());
    }

    @Test
    public void loadPeachJavaGenotypeFile() throws IOException
    {
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(PEACH_JAVA_GENOTYPES_FILE);

        assertEquals(3, peachGenotypes.size());

        assertEquals("GENE", peachGenotypes.get(0).gene());
        assertEquals("*1_HET", peachGenotypes.get(0).haplotype());
        assertEquals("*1", peachGenotypes.get(0).allele());
        assertEquals(1, peachGenotypes.get(0).alleleCount());
        assertEquals("Reduced Function", peachGenotypes.get(0).function());
        assertEquals("Capecitabine;otherDrug", peachGenotypes.get(0).linkedDrugs());
        assertEquals("link;otherLink", peachGenotypes.get(0).urlPrescriptionInfo());
        assertEquals("", peachGenotypes.get(0).panelVersion());
        assertEquals("", peachGenotypes.get(0).repoVersion());

        assertEquals("GENE2", peachGenotypes.get(1).gene());
        assertEquals("*3_HOM", peachGenotypes.get(1).haplotype());
        assertEquals("*3", peachGenotypes.get(1).allele());
        assertEquals(2, peachGenotypes.get(1).alleleCount());
        assertEquals("Normal Function", peachGenotypes.get(1).function());
        assertEquals("Capecitabine", peachGenotypes.get(1).linkedDrugs());
        assertEquals("link", peachGenotypes.get(1).urlPrescriptionInfo());
        assertEquals("", peachGenotypes.get(1).panelVersion());
        assertEquals("", peachGenotypes.get(1).repoVersion());

        assertEquals("GENE3", peachGenotypes.get(2).gene());
        assertEquals("Unresolved Haplotype", peachGenotypes.get(2).haplotype());
        assertEquals("Unresolved Haplotype", peachGenotypes.get(2).allele());
        assertEquals(2, peachGenotypes.get(2).alleleCount());
        assertEquals("Unknown Function", peachGenotypes.get(2).function());
        assertEquals("Capecitabine", peachGenotypes.get(2).linkedDrugs());
        assertEquals("link", peachGenotypes.get(2).urlPrescriptionInfo());
        assertEquals("", peachGenotypes.get(2).panelVersion());
        assertEquals("", peachGenotypes.get(2).repoVersion());
    }

    @Test
    public void createEmptyPeachJavaGenotypeOutputLines()
    {
        List<PeachGenotype> peachGenotypes = Collections.emptyList();

        List<String> expectedLines = List.of("gene\thaplotype\tcount\tfunction\tlinkedDrugs\tprescriptionUrls");

        assertEquals(expectedLines, PeachGenotypeFile.toLines(peachGenotypes));
    }

    @Test
    public void createNonEmptyPeachJavaGenotypeOutputLines()
    {
        List<PeachGenotype> peachGenotypes = List.of(
                ImmutablePeachGenotype.builder()
                        .gene("GENE")
                        .allele("*1")
                        .alleleCount(1)
                        .function("Reduced Function")
                        .linkedDrugs("Capecitabine;otherDrug")
                        .urlPrescriptionInfo("link;otherLink")
                        .panelVersion("V1")
                        .repoVersion("V2")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("GENE2")
                        .allele("*3")
                        .alleleCount(2)
                        .function("Normal Function")
                        .linkedDrugs("Capecitabine")
                        .urlPrescriptionInfo("link")
                        .panelVersion("")
                        .repoVersion("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("GENE3")
                        .allele("Unresolved Haplotype")
                        .alleleCount(2)
                        .function("Unknown Function")
                        .linkedDrugs("Capecitabine")
                        .urlPrescriptionInfo("link")
                        .build()
        );

        List<String> expectedLines = List.of(
                "gene\thaplotype\tcount\tfunction\tlinkedDrugs\tprescriptionUrls",
                "GENE\t*1\t1\tReduced Function\tCapecitabine;otherDrug\tlink;otherLink",
                "GENE2\t*3\t2\tNormal Function\tCapecitabine\tlink",
                "GENE3\tUnresolved Haplotype\t2\tUnknown Function\tCapecitabine\tlink"
        );

        assertEquals(expectedLines, PeachGenotypeFile.toLines(peachGenotypes));
    }
}
