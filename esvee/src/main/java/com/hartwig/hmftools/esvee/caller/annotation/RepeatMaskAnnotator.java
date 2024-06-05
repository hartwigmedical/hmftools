package com.hartwig.hmftools.esvee.caller.annotation;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarAlignedLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.gripss.RepeatMaskData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.esvee.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.caller.Variant;

import htsjdk.samtools.Cigar;

public class RepeatMaskAnnotator
{
    private final RepeatMaskAnnotations mAnnotationCache;

    private static final int MIN_OVERLAP = 20;
    private static final double MIN_COVERAGE_PERC = 0.1;
    private static final double POLY_A_T_PERC = 0.9;

    private static final RepeatMaskData POLY_T_DATA = new RepeatMaskData(
            0, new BaseRegion(0, 1), 0, ' ', "T(n)", "Simple_repeat");

    private static final RepeatMaskData POLY_A_DATA = new RepeatMaskData(
            0, new BaseRegion(0, 1), 0, ' ', "A(n)", "Simple_repeat");

    public RepeatMaskAnnotator()
    {
        mAnnotationCache = new RepeatMaskAnnotations();
    }

    public boolean hasData() { return mAnnotationCache.hasData(); }

    public void annotateVariants(final List<Variant> variantList)
    {
        int annotated = 0;
        for(Variant var : variantList)
        {
            if(var.insertSequence().isEmpty())
                continue;

            final String alignments = var.breakendStart().Context.getAttributeAsString(INSALN, "");
            if(alignments.isEmpty())
                continue;

            RepeatMaskAnnotation rmAnnotation = annotate(var.insertSequence(), alignments);

            if(rmAnnotation != null)
            {
                var.setRepeatMaskAnnotation(rmAnnotation);
                ++annotated;
            }
        }

        SV_LOGGER.debug("marked {} repeat mask annotations", annotated);
    }

    public RepeatMaskAnnotation annotate(final String insertSequence, final String alignmentsStr)
    {
        // List<AlignmentData> alignments = fromInsertSequenceAlignments(alignmentsStr);
        List<AlternativeAlignment> alignments = AlternativeAlignment.fromVcfTag(alignmentsStr);

        if(alignments == null || alignments.isEmpty())
            return null;

        double insSeqLength = insertSequence.length();

        // first check for a poly-A or T
        for(AlternativeAlignment alignment : alignments)
        {
            String matchedBases = extractMatchedBases(insertSequence, alignment);
            double polyAtPerc = calcPolyATPercent(matchedBases);

            if(polyAtPerc >= POLY_A_T_PERC)
            {
                double coverage = matchedBases.length() / insSeqLength;
                long aCount = matchedBases.chars().filter(x -> x == 'A').count();
                RepeatMaskData rmData = aCount > matchedBases.length() / 2 ? POLY_A_DATA : POLY_T_DATA;
                return new RepeatMaskAnnotation(rmData, coverage, alignment);
            }
        }

        RepeatMaskData topRm = null;
        AlternativeAlignment topAlignment = null;
        double topCoverage = 0;

        for(AlternativeAlignment alignment : alignments)
        {
            int cigarLength = cigarAlignedLength(alignment.cigar());
            int alignmentEnd = alignment.Position + cigarLength - 1;
            BaseRegion alignmentRegion = new BaseRegion(alignment.Position, alignmentEnd);
            List<RepeatMaskData> rmMatches = mAnnotationCache.findMatches(alignment.Chromosome, alignmentRegion);

            if(rmMatches.isEmpty())
                continue;

            for(RepeatMaskData rmData : rmMatches)
            {
                int overlap = min(alignmentRegion.end(), rmData.Region.end()) - max(alignmentRegion.start(), rmData.Region.start()) + 1;

                if(overlap < MIN_OVERLAP)
                    continue;

                double coverage = overlap / insSeqLength;

                if(coverage < MIN_COVERAGE_PERC)
                    continue;

                if(coverage > topCoverage)
                {
                    topCoverage = coverage;
                    topRm = rmData;
                    topAlignment = alignment;
                }
            }
        }

        return topRm != null ? new RepeatMaskAnnotation(topRm, topCoverage, topAlignment) : null;
    }

    private static String extractMatchedBases(final String insertSequence, final AlternativeAlignment alignment)
    {
        Cigar cigar = alignment.cigar();

        int matchStartPos = leftSoftClipped(cigar) ? cigar.getFirstCigarElement().getLength() : 0;
        int matchBases = cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() == M || x.getOperator() == I)
                .mapToInt(x -> x.getLength()).sum();

        int insSeqLength = insertSequence.length();
        int startPos = min(matchStartPos, insSeqLength);
        int endPos = min(matchStartPos + matchBases, insSeqLength);

        return insertSequence.substring(startPos, endPos);
    }

    private static double calcPolyATPercent(final String sequence)
    {
        int count = 0;

        for(int i = 0; i < sequence.length(); ++i)
        {
            if(sequence.charAt(i) == 'T' || sequence.charAt(i) == 'A')
                ++count;
        }

        return count / (double)sequence.length();
    }

    public boolean load(final String filename, final RefGenomeVersion refGenomeVersion)
    {
        return mAnnotationCache.load(filename, refGenomeVersion);
    }
}
