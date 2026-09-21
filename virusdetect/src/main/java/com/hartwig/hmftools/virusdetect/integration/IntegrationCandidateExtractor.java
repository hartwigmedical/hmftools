package com.hartwig.hmftools.virusdetect.integration;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_PAIRED_INSERT_LENGTH;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_SINGLE_INSERT_LENGTH;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.virusdetect.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

// Reads the ESVEE unfiltered VCF and keeps every SV which could be a viral integration.
// Uses the ESVEE unfiltered VCF because the viral integrations are interesting even if ESVEE decided to filter.
public class IntegrationCandidateExtractor
{
    private final String mTumorSampleId;

    private static final int NO_GENOTYPE_ORDINAL = -1;

    private static final Logger LOGGER = LogManager.getLogger(IntegrationCandidateExtractor.class);

    public IntegrationCandidateExtractor(String tumorSampleId)
    {
        mTumorSampleId = tumorSampleId;
    }

    public List<IntegrationCandidate> extract(String vcfFile)
    {
        try(VcfFileReader reader = new VcfFileReader(vcfFile))
        {
            if(!reader.fileValid())
            {
                throw new UserInputError("ESVEE VCF could not be read: " + vcfFile);
            }

            StructuralVariantFactory svFactory = StructuralVariantFactory.build(new AlwaysPassFilter());
            setGenotypeOrdinals(svFactory, reader, vcfFile);

            int variantCount = 0;
            for(VariantContext context : reader.iterator())
            {
                ++variantCount;
                svFactory.addVariantContext(context);
            }

            List<IntegrationCandidate> candidates = svFactory.results().stream()
                    .map(IntegrationCandidateExtractor::toCandidate)
                    .filter(Objects::nonNull)
                    .toList();

            LOGGER.info(
                    "Read {} variant records, {} breakends never paired, {} integration candidates",
                    variantCount, svFactory.unmatched().size(), candidates.size());

            return candidates;
        }
    }

    @Nullable
    private static IntegrationCandidate toCandidate(StructuralVariant variant)
    {
        StructuralVariantLeg endLeg = variant.end();
        String insertSequence = variant.insertSequence();
        
        int minInsertLength = endLeg == null ? MIN_SINGLE_INSERT_LENGTH : MIN_PAIRED_INSERT_LENGTH;
        if(insertSequence.length() < minInsertLength)
        {
            return null;
        }

        VariantContext startContext = variant.startContext();
        if(startContext == null)
        {
            throw new IllegalStateException("SV has no variant context: " + variant.id());
        }

        return new IntegrationCandidate(
                variant.id(),
                variant.type(),
                requireNonNull(variant.filter()),
                HostBreakend.from(variant.start()),
                endLeg != null ? HostBreakend.from(endLeg) : null,
                insertSequence,
                startContext.hasAttribute(LINE_SITE),
                InsertRepeat.from(variant),
                requireNonNull(variant.insertSequenceAlignments()));
    }

    // Fragment counts are attributed per sample, so the tumor genotype is resolved by name rather than by assuming an
    // ordinal. A second sample, if there is one, is the reference.
    private void setGenotypeOrdinals(StructuralVariantFactory svFactory, VcfFileReader reader, String vcfFile)
    {
        Map<String, Integer> ordinals = reader.genotypeOrdinals();
        Integer tumorOrdinal = ordinals.get(mTumorSampleId);
        if(tumorOrdinal == null)
        {
            throw new UserInputError(String.format(
                    "ESVEE VCF has no genotype for sample %s, found %s: %s", mTumorSampleId, ordinals.keySet(), vcfFile));
        }

        int referenceOrdinal = ordinals.size() == 2
                ? ordinals.values().stream().filter(ordinal -> ordinal != tumorOrdinal.intValue()).findFirst().orElseThrow()
                : NO_GENOTYPE_ORDINAL;

        svFactory.setGenotypeOrdinals(referenceOrdinal, tumorOrdinal);
    }
}
