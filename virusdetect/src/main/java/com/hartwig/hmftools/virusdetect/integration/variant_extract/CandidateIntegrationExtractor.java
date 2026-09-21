package com.hartwig.hmftools.virusdetect.integration.variant_extract;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.INTEGRATION_SGL_INSERT_LENGTH_MIN;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.INTEGRATION_VARIANT_INSERT_LENGTH_MIN;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.virusdetect.common.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

// Reads the ESVEE unfiltered VCF and keeps every SV which could be a viral integration.
// Uses the ESVEE unfiltered VCF because the viral integrations are interesting even if ESVEE decided to filter.
public class CandidateIntegrationExtractor
{
    private final String mTumorSampleId;

    private static final int NO_GENOTYPE_ORDINAL = -1;

    private static final Logger LOGGER = LogManager.getLogger(CandidateIntegrationExtractor.class);

    public CandidateIntegrationExtractor(String tumorSampleId)
    {
        mTumorSampleId = tumorSampleId;
    }

    public List<CandidateIntegration> extract(String vcfFile)
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

            List<CandidateIntegration> candidates = svFactory.results().stream()
                    .map(CandidateIntegrationExtractor::toCandidate)
                    .filter(Objects::nonNull)
                    .toList();

            LOGGER.info(
                    "Read {} variant records, {} breakends never paired, {} integration candidates",
                    variantCount, svFactory.unmatched().size(), candidates.size());

            return candidates;
        }
    }

    @Nullable
    private static CandidateIntegration toCandidate(StructuralVariant variant)
    {
        StructuralVariantLeg endLeg = variant.end();
        String insertSequence = variant.insertSequence();

        int minInsertLength = endLeg == null ? INTEGRATION_SGL_INSERT_LENGTH_MIN : INTEGRATION_VARIANT_INSERT_LENGTH_MIN;
        if(insertSequence.length() < minInsertLength)
        {
            return null;
        }

        VariantContext startContext = variant.startContext();
        if(startContext == null)
        {
            throw new IllegalStateException("SV has no variant context: " + variant.id());
        }

        return new CandidateIntegration(
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
