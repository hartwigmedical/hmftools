package com.hartwig.hmftools.virusdetect.integration;

import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_PAIRED_INSERT_LENGTH;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_SINGLE_INSERT_LENGTH;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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

// Reads the ESVEE unfiltered VCF and keeps every SV carrying an inserted sequence long enough to be worth aligning to
// the viral reference. The unfiltered VCF is read because viral junctions are often single-sided and low support,
// which is what ESVEE's filters remove; the filter string is carried through as data instead.
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

            // The factory buffers a paired variant until both its mate records arrive, and drops non-human contigs, which
            // viral sequence never reaches anyway: ESVEE reports it only as inserted sequence in the ALT allele.
            StructuralVariantFactory svFactory = StructuralVariantFactory.build(new AlwaysPassFilter());
            setGenotypeOrdinals(svFactory, reader, vcfFile);

            List<IntegrationCandidate> candidates = new ArrayList<>();
            int variantCount = 0;

            for(VariantContext context : reader.iterator())
            {
                ++variantCount;

                int completedBefore = svFactory.results().size();
                svFactory.addVariantContext(context);
                if(svFactory.results().size() == completedBefore)
                {
                    continue;
                }

                StructuralVariant variant = svFactory.results().remove(completedBefore);
                IntegrationCandidate candidate = toCandidate(variant);
                if(candidate != null)
                {
                    candidates.add(candidate);
                }
            }

            LOGGER.info(
                    "read {} variant records, {} breakends never paired, {} integration candidates",
                    variantCount, svFactory.unmatched().size(), candidates.size());

            return candidates;
        }
    }

    // Insert length is the only gate in this phase. Everything past it reaches the output, aligned to a virus or not.
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
                variant.filter() != null ? variant.filter() : "",
                HostBreakend.from(variant.start()),
                endLeg != null ? HostBreakend.from(endLeg) : null,
                insertSequence,
                startContext.hasAttribute(LINE_SITE),
                variant.insertSequenceRepeatClass(),
                variant.insertSequenceRepeatType(),
                variant.insertSequenceRepeatOrientation(),
                variant.insertSequenceRepeatCoverage(),
                variant.insertSequenceAlignments());
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
