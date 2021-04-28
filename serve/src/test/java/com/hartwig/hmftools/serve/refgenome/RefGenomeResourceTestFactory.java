package com.hartwig.hmftools.serve.refgenome;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public final class RefGenomeResourceTestFactory {

    private static final String REF_GENOME_37_FASTA_FILE = Resources.getResource("refgenome/v37/ref.fasta").getPath();

    private RefGenomeResourceTestFactory() {
    }

    @NotNull
    public static RefGenomeResource buildTest37() {
        IndexedFastaSequenceFile refSequence;
        try {
            refSequence = new IndexedFastaSequenceFile(new File(REF_GENOME_37_FASTA_FILE));
        } catch (FileNotFoundException e) {
            throw new IllegalStateException("Could not create ref sequence from " + REF_GENOME_37_FASTA_FILE);
        }

        return ImmutableRefGenomeResource.builder()
                .refSequence(refSequence)
                .driverGenes(Lists.newArrayList())
                .knownFusionCache(new KnownFusionCache())
                .canonicalTranscriptPerGeneMap(HmfGenePanelSupplier.allGenesMap37())
                .proteinResolver(new TestProteinResolver())
                .build();
    }

    private static class TestProteinResolver implements ProteinResolver {

        @NotNull
        @Override
        public List<VariantHotspot> resolve(@NotNull final String gene, @Nullable final String specificTranscript,
                @NotNull final String proteinAnnotation) {
            return Lists.newArrayList(ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build());
        }

        @NotNull
        @Override
        public Set<String> unresolvedProteinAnnotations() {
            return Sets.newHashSet();
        }
    }
}
