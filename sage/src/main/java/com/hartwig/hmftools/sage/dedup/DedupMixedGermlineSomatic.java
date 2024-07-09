package com.hartwig.hmftools.sage.dedup;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.filter.SoftFilter.DEDUP_MIXED_GERMLINE_SOMATIC;
import static com.hartwig.hmftools.sage.filter.SoftFilterConfig.getTieredSoftFilterConfig;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.filter.SoftFilter;
import com.hartwig.hmftools.sage.filter.SoftFilterConfig;
import com.hartwig.hmftools.sage.select.TranscriptSelector;

public class DedupMixedGermlineSomatic
{
    private final TranscriptSelector mSelector;
    private final FilterConfig mFilterConfig;
    private int mNextMixedId;

    private static int MAX_DISTANCE = 5;

    public DedupMixedGermlineSomatic(final List<TranscriptData> transcripts, final FilterConfig filterConfig)
    {
        mSelector = new TranscriptSelector(transcripts);
        mFilterConfig = filterConfig;
        mNextMixedId = 1;
    }

    public void dedupVariants(final List<SageVariant> variants)
    {
        /*
        previous logic

        3. MixedSomaticGermlineIdentifier:
        - only applies to MNVs and SNVs
        - max distance of 10 bases
        - only applies to passing MNVs and SNVs which are germline soft-filtered (ie no passing, and not tumor filtered)
        - looks for a combination of the above 2 conditions
        - if the MNV contains the SNV, then record this against both variants

        4. MixedSomaticGermlineDedup
        - max distance of 10 bases
        - looks for passing SNV and mixed germline MNV from step 3 above
        - must have a common LPS
        - if the MNV contains the SNV and they fall in the same codon (this was a chance of being wrong for -ve strand as I've mentioned since assumed codons started on transcript start) then
        - mark the MNV as dedup if has more differences within its bases that overlap the codon, otherwise mark the SNV as dedup

        */

        List<SageVariant> candidates = variants.stream()
                .filter(x -> isGermlineFilteredSnv(x) || isPassingSnv(x) || isPassingMnv(x))
                .collect(Collectors.toList());

        int index = 0;
        while(index < candidates.size() - 1)
        {
            SageVariant firstVariant = candidates.get(index);

            int maxPos = firstVariant.position() + firstVariant.alt().length() - 1;
            int nextIndex = index + 1;

            SageVariant mnv = null;
            SageVariant snvPassing = null;
            SageVariant snvGermline = null;

            if(isPassingMnv(firstVariant))
                mnv = firstVariant;
            else if(isPassingSnv(firstVariant))
                snvPassing = firstVariant;
            else
                snvGermline = firstVariant;

            while(nextIndex < candidates.size())
            {
                SageVariant nextVariant = candidates.get(nextIndex);

                if(nextVariant.position() - maxPos > MAX_DISTANCE)
                    break;

                if(mnv == null && isPassingMnv(nextVariant))
                    mnv = nextVariant;
                else if(snvPassing == null && isPassingSnv(nextVariant))
                    snvPassing = nextVariant;
                else if(snvGermline == null && isGermlineFilteredSnv(nextVariant))
                    snvGermline = nextVariant;

                // check the combo
                if(mnv != null && snvPassing != null && snvGermline != null)
                {
                    checkGroup(mnv, snvPassing, snvGermline);
                    ++nextIndex; // skip past the last in the group
                    break;
                }

                ++nextIndex;
            }

            ++index;
        }
    }

    private void checkGroup(final SageVariant mnv, final SageVariant snvPassing, final SageVariant snvGermline)
    {
        if(!snvPassing.hasMatchingLps(mnv.localPhaseSets()))
            return;

        if(VariantDeduper.longerContainsShorter(snvGermline, mnv))
        {
            if(mnv.mixedGermlineImpact() == 0)
            {
                mnv.mixedGermlineImpact(mNextMixedId++);
            }

            snvGermline.mixedGermlineImpact(mnv.mixedGermlineImpact());

            if(VariantDeduper.longerContainsShorter(snvPassing, mnv))
            {
                final BaseRegion positionCodon = findCodon(snvPassing.position());

                if(positionCodon != null && keepMnv(positionCodon, snvPassing.variant(), mnv.variant()))
                {
                    snvPassing.filters().add(DEDUP_MIXED_GERMLINE_SOMATIC);
                }
                else
                {
                    mnv.filters().add(DEDUP_MIXED_GERMLINE_SOMATIC);
                }
            }
        }
    }

    public static boolean keepMnv(final BaseRegion codon, final SimpleVariant somaticSnv, final SimpleVariant mixedMnv)
    {
        int snvCodonDifferences = codonDifferences(codon, somaticSnv);
        int mnvCodonDifferences = codonDifferences(codon, mixedMnv);

        return mnvCodonDifferences > snvCodonDifferences;
    }

    private static boolean isGermlineFilteredSnv(final SageVariant variant)
    {
        return variant.isSnv() && !variant.isPassing() && SoftFilter.isGermlineAndNotTumorFiltered(variant.filters());
    }

    private static boolean isPassingMnv(final SageVariant variant)
    {
        return variant.isMnv() && variant.isPassing() && variant.hasLocalPhaseSets();
    }

    private static boolean isPassingSnv(final SageVariant variant)
    {
        return variant.isSnv() && variant.isPassing() && variant.hasLocalPhaseSets();
    }

    public static int codonDifferences(final BaseRegion codon, final SimpleVariant variant)
    {
        if(codon.start() > variant.end() || variant.Position > codon.end())
            return 0;

        int overlapStart = max(codon.start(), variant.Position);
        int overlapEnd = min(codon.end(), variant.end());

        int difference = 0;
        for(int position = overlapStart; position <= overlapEnd; position++)
        {
            int variantIndex = position - variant.position();

            if(variant.Ref.charAt(variantIndex) != variant.Alt.charAt(variantIndex))
            {
                difference++;
            }
        }
        return difference;
    }

    private BaseRegion findCodon(int position)
    {
        final TranscriptData transcript = mSelector.select(position);

        if(transcript == null || transcript.nonCoding())
            return null;

        if(!positionWithin(position, transcript.CodingStart, transcript.CodingEnd))
            return null;

        // find the codon surrounding this position
        ExonData exon = transcript.exons().stream().filter(x -> positionWithin(position, x.Start, x.End)).findFirst().orElse(null);

        if(exon == null)
            return null;

        int codingPhase = calcExonicCodingPhase(exon, transcript.CodingStart, transcript.CodingEnd, transcript.Strand, position);
        int codonStart = position;
        int codonEnd = position;

        if(transcript.Strand == POS_STRAND)
        {
            if(codingPhase == PHASE_1)
            {
                codonEnd = min(exon.End, position + 2);
            }
            else if(codingPhase == PHASE_2)
            {
                codonStart = max(exon.Start, position - 1);
                codonEnd = min(exon.End, position + 1);
            }
            else if(codingPhase == PHASE_0)
            {
                codonStart = max(exon.Start, position - 2);
            }
        }
        else
        {
            if(codingPhase == PHASE_1)
            {
                codonStart = max(exon.Start, position - 2);
            }
            else if(codingPhase == PHASE_2)
            {
                codonStart = max(exon.Start, position - 1);
                codonEnd = min(exon.End, position + 1);
            }
            else if(codingPhase == PHASE_0)
            {
                codonEnd = min(exon.End, position + 2);
            }
        }

        return new BaseRegion(codonStart, codonEnd);
    }

}
