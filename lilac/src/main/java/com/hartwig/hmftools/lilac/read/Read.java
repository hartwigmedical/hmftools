package com.hartwig.hmftools.lilac.read;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.generateMappedCoords;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractLowQualIndices;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.lilac.LilacConfig.isSbx;
import static com.hartwig.hmftools.lilac.LilacConfig.isUltima;
import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_TRIM_PERC;
import static com.hartwig.hmftools.lilac.LilacUtils.belowMinQual;
import static com.hartwig.hmftools.lilac.ReferenceData.refBases;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class Read
{
    public final String Id;
    public final int SoftClippedStart; // soft-clipped bases at start and end
    public final int SoftClippedEnd;

    // coords which cover the coding region of the applicable gene
    public final int PositionStart; // adjusted by soft-clipped bases
    public final int PositionEnd;

    // corresponding read base indices
    public final int ReadIndexStart;
    public final int ReadIndexEnd;

    private final List<Indel> mValidIndels; // non-ignored indels within the coding region
    private final List<Indel> mIgnoredIndels; // ignored indels within the coding region
    private final SAMRecord mRecord;
    private final int mTrimmedBases;

    public Read(
            final String id, int softClippedStart, int softClippedEnd, final List<Indel> indels, int positionStart,
            int positionEnd, int readIndexStart, int readIndexEnd, final SAMRecord record, int trimmedBases)
    {
        Id = id;
        SoftClippedStart = softClippedStart;
        SoftClippedEnd = softClippedEnd;
        mValidIndels = indels.stream().filter(x -> !x.Ignore).toList();
        mIgnoredIndels = indels.stream().filter(x -> x.Ignore).toList();
        PositionStart = positionStart;
        PositionEnd = positionEnd;
        ReadIndexStart = readIndexStart;
        ReadIndexEnd = readIndexEnd;
        mRecord = record;
        mTrimmedBases = trimmedBases;
    }

    public String readInfo()
    {
        return format("%s:%d-%d %s",
                mRecord.getContig(), mRecord.getStart(), mRecord.getEnd(), mRecord.getCigarString());
    }

    public String toString()
    {
        return format("%s %s indels(%d) ignoredIndels(%d) sc(%d/%d)",
                Id, readInfo(), mValidIndels.size(), mIgnoredIndels.size(), SoftClippedStart, SoftClippedEnd);
    }

    public List<Indel> getValidIndels() { return mValidIndels; }
    public List<Indel> getIgnoredIndels() { return mIgnoredIndels; }
    public int trimmedBases() { return mTrimmedBases; }

    public SAMRecord bamRecord() { return mRecord; }

    public int maxValidIndelSize()
    {
        return mValidIndels.stream().mapToInt(x -> abs(x.Length)).max().orElse(0);
    }

    public final boolean containsSoftClip() { return SoftClippedStart > 0 || SoftClippedEnd > 0; }

    public boolean containsValidIndel() { return !mValidIndels.isEmpty(); }

    public void populateCodingRegion(final char[] readBases, final byte[] readQuals, boolean reverseCompliment)
    {
        Set<Integer> ultimaLowQualIndices = isUltima() ? Sets.newHashSet(extractLowQualIndices(mRecord)) : null;
        if(reverseCompliment)
        {
            int index = readBases.length - 1;
            for(int i = ReadIndexStart; i <= ReadIndexEnd; ++i)
            {
                readBases[index] = Nucleotides.swapDnaBase(mRecord.getReadString().charAt(i));
                byte qual;
                if(ultimaLowQualIndices == null)
                    qual = mRecord.getBaseQualities()[i];
                else
                    qual = ultimaLowQualIndices.contains(i) ? (byte) (LOW_BASE_QUAL_THRESHOLD - 1) : LOW_BASE_QUAL_THRESHOLD;


                readQuals[index] = qual;
                --index;
            }

        }
        else
        {
            int index = 0;
            for(int i = ReadIndexStart; i <= ReadIndexEnd; ++i)
            {
                readBases[index] = mRecord.getReadString().charAt(i);
                byte qual;
                if(ultimaLowQualIndices == null)
                    qual = mRecord.getBaseQualities()[i];
                else
                    qual = ultimaLowQualIndices.contains(i) ? (byte) (LOW_BASE_QUAL_THRESHOLD - 1) : LOW_BASE_QUAL_THRESHOLD;

                readQuals[index] = qual;
                ++index;
            }
        }
    }

    public List<Read> alignmentsOnly()
    {
        final String chromosome = mRecord.getContig();
        final ChrBaseRegion outerRegion = new ChrBaseRegion(chromosome, PositionStart, PositionEnd);

        return mRecord.getAlignmentBlocks().stream()
                .map(x -> new ChrBaseRegion(chromosome, x.getReferenceStart(), x.getReferenceStart() + x.getLength()))
                .filter(x -> outerRegion.overlaps(x))
                .map(x -> new BaseRegion(max(outerRegion.start(), x.start()), min(outerRegion.end(), x.end())))
                .map(x -> createRead(x, mRecord, false, false))
                .collect(Collectors.toList());
    }

    public static Read createRead(final BaseRegion codingRegion, final SAMRecord record, boolean includeSoftClips, boolean includeIndels)
    {
        int scLengthLeft = leftSoftClipLength(record);
        int scLengthRight = rightSoftClipLength(record);
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        // initialise positions and read indices to non-soft-clipped alignments
        int positionStart = alignmentStart;
        int positionEnd = alignmentEnd;

        int readIndexStart = scLengthLeft;
        int readIndexEnd = record.getReadBases().length - scLengthRight - 1;

        int softClipStart = scLengthLeft;
        int softClipEnd = scLengthRight;

        // ignore soft-clip positions on the 3' end which run past the mate's 5' end for overlapping fragments
        if(record.getReadPairedFlag() && abs(record.getInferredInsertSize()) <= record.getReadBases().length)
        {
            if(record.getReadNegativeStrandFlag())
                softClipStart = 0;
            else
                softClipEnd = 0;
        }

        int lowQualTrimmedBases = findLowQualTrimmedBases(record);

        if(lowQualTrimmedBases > 0)
        {
            // trim or cancel any soft-clip bases from low BQ trimming
            if(record.getReadNegativeStrandFlag())
            {
                softClipStart = max(softClipStart - lowQualTrimmedBases, 0);

                readIndexStart += (lowQualTrimmedBases - scLengthLeft);
                positionStart = record.getReferencePositionAtReadPosition(readIndexStart + 1);
            }
            else
            {
                softClipEnd = max(softClipEnd - lowQualTrimmedBases, 0);

                readIndexEnd -= (lowQualTrimmedBases - scLengthRight);
                positionEnd = record.getReferencePositionAtReadPosition(readIndexEnd + 1);
            }
        }

        int unclippedStart = softClipStart > 0 ? max(alignmentStart - softClipStart, codingRegion.start()) : 0;
        int unclippedEnd = softClipEnd > 0 ? min(alignmentEnd + softClipEnd, codingRegion.end()) : 0;

        // include soft-clip bases if permitted
        if(includeSoftClips)
        {
            if(softClipStart > 0)
                positionStart = unclippedStart;

            if(softClipEnd > 0)
                positionEnd = unclippedEnd;
        }

        // restrict to coding region
        positionStart = max(positionStart, codingRegion.start());
        positionEnd = min(positionEnd, codingRegion.end());

        // use our bespoke method since it allows positions in the soft-clipped region
        readIndexStart = getReadIndexFromPosition(
                alignmentStart, record.getCigar().getCigarElements(), positionStart, false, true);

        readIndexEnd = getReadIndexFromPosition(
                alignmentStart, record.getCigar().getCigarElements(), positionEnd, false, true);

        if(record.getCigar().containsOperator(CigarOperator.N))
        {
            // handle splits which do not align with exon boundaries (which the coding region represents)
            List<int[]> mappedCoords = generateMappedCoords(record.getCigar(), record.getAlignmentStart());

            for(int[] mappedCoord : mappedCoords)
            {
                if(mappedCoord[SE_END] < codingRegion.start())
                    continue;

                if(mappedCoord[SE_START] > codingRegion.end())
                    break;

                positionStart = max(max(positionStart, mappedCoord[SE_START]), codingRegion.start());
                positionEnd = min(min(positionEnd, mappedCoord[SE_END]), codingRegion.end());

                readIndexStart = record.getReadPositionAtReferencePosition(positionStart, true) - 1;
                readIndexEnd = record.getReadPositionAtReferencePosition(positionEnd, true) - 1;
                break;
            }
        }

        int softClippedStart = max(alignmentStart - positionStart, 0);
        int softClippedEnd = max(0, positionEnd - alignmentEnd);

        List<Indel> indels = includeIndels ? indels(positionStart, positionEnd, record) : Lists.newArrayList();

        return new Read(
                record.getReadName(), softClippedStart, softClippedEnd, indels,
                positionStart, positionEnd, readIndexStart, readIndexEnd, record, lowQualTrimmedBases);
    }

    private static List<Indel> indels(int startPosition, int endPosition, final SAMRecord record)
    {
        final List<Indel> indels = Lists.newArrayList();
        final Set<Integer> ultimaLowQualIndices = isUltima() ? Sets.newHashSet(extractLowQualIndices(record)) : null;

        CigarHandler cigarHandler = new CigarHandler()
        {
            public void handleInsert(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                if(startPosition <= refPosition && refPosition <= endPosition)
                {
                    char baseChar = record.getReadString().charAt(readIndex);
                    boolean ignore = false;
                    if(isUltima() || isSbx())
                    {
                        for(int i = readIndex + 1; i <= readIndex + element.getLength(); i++)
                        {
                            boolean isLowQual = ultimaLowQualIndices != null
                                    ? ultimaLowQualIndices.contains(i)
                                    : belowMinQual(record.getBaseQualities()[i]);
                            if(isLowQual)
                            {
                                ignore = true;
                                break;
                            }
                        }
                    }

                    indels.add(new Indel(record.getContig(), refPosition, readIndex, String.valueOf(baseChar),
                            record.getReadString().substring(readIndex, readIndex + element.getLength() + 1), ignore));
                }
            }

            public void handleDelete(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                if(startPosition <= refPosition && refPosition <= endPosition)
                {
                    char baseChar = record.getReadString().charAt(readIndex);
                    int delPosStart = refPosition + 1;
                    int delPosEnd = delPosStart + element.getLength() - 1;
                    String delBases = baseChar + refBases(delPosStart, delPosEnd);
                    boolean ignore = false;
                    if(isUltima() || isSbx())
                    {
                        int leftIndex = readIndex;
                        boolean isLeftLowQual = ultimaLowQualIndices != null
                                ? ultimaLowQualIndices.contains(leftIndex)
                                : belowMinQual(record.getBaseQualities()[leftIndex]);
                        int rightIndex = leftIndex + 1;
                        boolean isRightLowQual = ultimaLowQualIndices != null
                                ? ultimaLowQualIndices.contains(rightIndex)
                                : belowMinQual(record.getBaseQualities()[rightIndex]);
                        if(isLeftLowQual || isRightLowQual)
                            ignore = true;
                    }

                    indels.add(new Indel(record.getContig(), refPosition, readIndex, delBases, String.valueOf(baseChar), ignore));
                }
            }
        };

        CigarHandler.traverseCigar(record, cigarHandler);
        return indels;
    }

    private static int findLowQualTrimmedBases(final SAMRecord record)
    {
        boolean fromStart = record.getReadNegativeStrandFlag();
        int baseIndex = fromStart ? 0 : record.getReadBases().length - 1;

        double lowestScore = 0;
        double currentScore = 0;
        int lastLowestScoreIndex = -1;

        Set<Integer> ultimaLowQualIndices = isUltima() ? Sets.newHashSet(extractLowQualIndices(record)) : null;
        while(baseIndex >= 0 && baseIndex < record.getReadBases().length)
        {
            boolean isLowQual = ultimaLowQualIndices != null ? ultimaLowQualIndices.contains(baseIndex) : belowMinQual(record.getBaseQualities()[baseIndex]);
            if(isLowQual)
            {
                currentScore -= LOW_QUAL_SCORE;

                if(currentScore <= lowestScore)
                {
                    lowestScore = currentScore;
                    lastLowestScoreIndex = baseIndex;
                }
            }
            else
            {
                ++currentScore;
            }

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lastLowestScoreIndex < 0)
            return 0;

        return fromStart ? lastLowestScoreIndex + 1 : record.getReadBases().length - lastLowestScoreIndex;
    }

    protected static final double LOW_QUAL_SCORE = 1 / LOW_BASE_TRIM_PERC - 1;
}
