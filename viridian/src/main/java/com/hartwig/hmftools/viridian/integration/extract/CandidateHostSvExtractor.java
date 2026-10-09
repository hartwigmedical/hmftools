package com.hartwig.hmftools.viridian.integration.extract;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_INSERT_LENGTH_MIN;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

// Reads the ESVEE unfiltered VCF and keeps every SV which could be a viral integration.
// Uses the ESVEE unfiltered VCF because the viral integrations are interesting even if ESVEE decided to filter.
public class CandidateHostSvExtractor
{
    private final String mTumorSampleId;

    private static final int NO_GENOTYPE_ORDINAL = -1;

    private static final Logger LOGGER = LogManager.getLogger(CandidateHostSvExtractor.class);

    public CandidateHostSvExtractor(String tumorSampleId)
    {
        mTumorSampleId = tumorSampleId;
    }

    public List<CandidateHostSv> extract(String vcfFile)
    {
        try(VcfFileReader reader = new VcfFileReader(vcfFile))
        {
            if(!reader.fileValid())
            {
                throw new UserInputError("ESVEE VCF could not be read: " + vcfFile);
            }

            StructuralVariantFactory svFactory = StructuralVariantFactory.build(new AlwaysPassFilter());
            setGenotypeOrdinals(svFactory, reader, vcfFile);

            for(VariantContext context : reader.iterator())
            {
                svFactory.addVariantContext(context);
            }

            List<CandidateHostSv> candidates = svFactory.results().stream()
                    .map(CandidateHostSvExtractor::toCandidate)
                    .filter(Objects::nonNull)
                    .toList();

            int unpairedBreakends = svFactory.unmatched().size();
            if(unpairedBreakends > 0)
            {
                LOGGER.warn("{} breakends never paired", unpairedBreakends);
            }

            return candidates;
        }
    }

    @Nullable
    private static CandidateHostSv toCandidate(StructuralVariant variant)
    {
        StructuralVariantLeg endLeg = variant.end();
        String insertSequence = variant.insertSequence();

        if(insertSequence.length() < INTEGRATION_INSERT_LENGTH_MIN)
        {
            return null;
        }

        VariantContext startContext = variant.startContext();
        if(startContext == null)
        {
            throw new IllegalStateException("SV has no variant context: " + variant.id());
        }

        return new CandidateHostSv(
                variant.type(),
                requireNonNull(variant.filter()),
                HostBreakend.from(startContext.getID(), variant.start()),
                endLeg != null ? HostBreakend.from(requireNonNull(variant.endContext()).getID(), endLeg) : null,
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
