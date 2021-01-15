package com.hartwig.hmftools.serve.extraction.hotspot;

import static org.junit.Assert.assertNotNull;

import java.io.FileNotFoundException;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

import org.junit.Test;

public class ProteinResolverFactoryTest {

    private static final String REF_GENOME_FASTA_FILE = Resources.getResource("refgenome/ref.fasta").getPath();

    @Test
    public void canCreateTransvarResolver() throws FileNotFoundException {
        assertNotNull(ProteinResolverFactory.transvarWithRefGenome(RefGenomeVersion.V37, REF_GENOME_FASTA_FILE, Maps.newHashMap()));
    }

}