package com.hartwig.hmftools.serve.hotspot;

import java.io.FileNotFoundException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.transvar.Transvar;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ProteinResolverFactory {

    private static final Logger LOGGER = LogManager.getLogger(ProteinResolverFactory.class);

    private ProteinResolverFactory() {
    }

    @NotNull
    public static ProteinResolver transvarWithRefGenome(@NotNull RefGenomeVersion refGenomeVersion,
            @NotNull String refGenomeFastaFile) throws FileNotFoundException {
        LOGGER.info("Creating protein to hotspot resolver with ref genome version '{}' and fasta path '{}'",
                refGenomeVersion,
                refGenomeFastaFile);
        return Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile);
    }

    @NotNull
    public static ProteinResolver dummy() {
        return new ProteinResolver() {
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
        };
    }
}
