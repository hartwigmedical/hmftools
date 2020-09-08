package com.hartwig.hmftools.serve.hotspot;

import java.io.FileNotFoundException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.transvar.Transvar;
import com.hartwig.hmftools.vicc.util.DetermineHotspot;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ProteinToHotspotConverter {

    private static final Logger LOGGER = LogManager.getLogger(ProteinToHotspotConverter.class);

    private static final int MAX_INFRAME_BASE_LENGTH = 50;

    @NotNull
    private final ProteinResolver proteinResolver;

    @NotNull
    public static ProteinToHotspotConverter transvarWithRefGenome(@NotNull RefGenomeVersion refGenomeVersion,
            @NotNull String refGenomeFastaFile) throws FileNotFoundException {
        LOGGER.info("Creating protein to hotspot resolver with ref genome version '{}' and fasta path '{}'",
                refGenomeVersion,
                refGenomeFastaFile);
        return new ProteinToHotspotConverter(Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile));
    }

    @NotNull
    public static ProteinToHotspotConverter dummy() {
        return new ProteinToHotspotConverter(new ProteinResolver() {
            @NotNull
            @Override
            public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull final String gene,
                    @Nullable final String specificTranscript, @NotNull final String proteinAnnotation) {
                return Lists.newArrayList();
            }

            @NotNull
            @Override
            public Set<String> unresolvedProteinAnnotations() {
                return Sets.newHashSet();
            }
        });
    }

    private ProteinToHotspotConverter(@NotNull ProteinResolver proteinResolver) {
        this.proteinResolver = proteinResolver;
    }

    @NotNull
    public List<VariantHotspot> resolveProteinAnnotation(@NotNull String gene, @Nullable String specificTranscript,
            @NotNull String proteinAnnotation) {
        if (DetermineHotspot.isResolvableProteinAnnotation(proteinAnnotation)) {
            return proteinResolver.extractHotspotsFromProteinAnnotation(gene, specificTranscript, proteinAnnotation);
        }

        return Lists.newArrayList();
    }

    @NotNull
    public Set<String> unresolvedProteinAnnotations() {
        return proteinResolver.unresolvedProteinAnnotations();
    }
}
