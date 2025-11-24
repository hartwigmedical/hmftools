package com.hartwig.hmftools.esvee.assembly.alignment;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DECOY_MAX_MISMATCHES;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DECOY_MIN_SCORE_FACTOR;

import static htsjdk.samtools.CigarOperator.M;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.WriteType;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

import htsjdk.samtools.Cigar;

public class AlignmentChecker
{
    private final AssemblyConfig mConfig;
    private final BwaMemAligner mDecoyAligner;
    private final BwaMemAligner mAligner;
    private int mSequenceCount;
    private int mSequenceMatched;
    private final BufferedWriter mWriter;

    public AlignmentChecker(final AssemblyConfig config, final BufferedWriter writer)
    {
        mConfig = config;
        mSequenceCount = 0;
        mSequenceMatched = 0;
        mWriter = writer;

        mAligner = config.AssemblyMapQualThreshold > 0 ? new BwaMemAligner(new BwaMemIndex(mConfig.RefGenomeImageFile)) : null;

        mDecoyAligner = config.DecoyGenome != null ? new BwaMemAligner(new BwaMemIndex(mConfig.DecoyGenome)) : null;
    }

    public int sequenceCount() { return mSequenceCount; }
    public int sequenceMatched() { return mSequenceMatched; };

    public boolean matchesDecoy(final JunctionAssembly assembly)
    {
        if(mDecoyAligner == null)
            return false;

        ++mSequenceCount;

        String fullSequence = assembly.formFullSequence();
        List<BwaMemAlignment> alignmentResults = mDecoyAligner.alignSeqs(List.of(fullSequence.getBytes())).get(0);

        if(alignmentResults.isEmpty())
            return false;

        BwaMemAlignment topAlignment = alignmentResults.get(0);

        if(topAlignment.getAlignerScore() == 0 || topAlignment.getCigar().isEmpty())
            return false;

        BwaMemAlignment nextAlignment = alignmentResults.size() > 1 ? alignmentResults.get(1) : null;

        Cigar cigar = CigarUtils.cigarFromStr(topAlignment.getCigar());
        int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        writeMatchInfo(mWriter, assembly, topAlignment, alignedBases, nextAlignment);

        int sequenceLength = fullSequence.length();
        int unalignedBases = sequenceLength - alignedBases;
        double scoreFactor = topAlignment.getAlignerScore() / (double)sequenceLength;

        if(unalignedBases <= DECOY_MAX_MISMATCHES && scoreFactor >= DECOY_MIN_SCORE_FACTOR)
        {
            ++mSequenceMatched;
            return true;
        }

        return false;
    }

    public boolean failsMappability(final JunctionAssembly assembly)
    {
        if(mAligner == null)
            return false;

        if(mConfig.AssemblyMapQualThreshold == 0 || averageAssemblySupportMapQual(assembly) >= mConfig.AssemblyMapQualThreshold)
            return false;

        // realign the ref base sequence and exclude if the results are too varied
        String refBaseSequence = assembly.formRefBaseSequence();
        List<BwaMemAlignment> alignmentResults = mAligner.alignSeqs(List.of(refBaseSequence.getBytes())).get(0);

        if(alignmentResults.isEmpty())
            return true;

        for(BwaMemAlignment alignment : alignmentResults)
        {
            if(alignment.getMapQual() >= mConfig.AssemblyMapQualThreshold)
                return false;

            if(alignment.getXATag() != null)
                return false;
        }

        return true;
    }

    private double averageAssemblySupportMapQual(final JunctionAssembly assembly)
    {
        int totalMapQual = assembly.support().stream().mapToInt(x -> x.mapQual()).sum();
        return assembly.supportCount() > 0 ? totalMapQual / (double)assembly.supportCount() : 0;
    }

    public static BufferedWriter initialiseWriter(final AssemblyConfig config)
    {
        if(!config.WriteTypes.contains(WriteType.DECOY_MATCHES))
            return null;

        if(config.OutputDir == null)
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(config.outputFilename(WriteType.DECOY_MATCHES));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("AssemblyChromosome");
            sj.add("AssemblyPosition");
            sj.add("AssemblyOrientation");
            sj.add("SequenceLength");
            sj.add("AlignCigar");
            sj.add("AlignMatchedBases");
            sj.add("AlignMismatches");
            sj.add("AlignLocation");
            sj.add("AlignScore");
            sj.add("AlignNextScore");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise decoy match writer: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeMatchInfo(
            final BufferedWriter writer, final JunctionAssembly assembly, final BwaMemAlignment topAlignment,
            final int alignedBases, final BwaMemAlignment nextAlignment)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(assembly.junction().Chromosome);
            sj.add(String.valueOf(assembly.junction().Position));
            sj.add(String.valueOf(assembly.junction().Orient));
            sj.add(String.valueOf(assembly.extensionLength() + assembly.refBaseLength()));
            sj.add(topAlignment.getCigar());
            sj.add(String.valueOf(alignedBases));
            sj.add(String.valueOf(topAlignment.getNMismatches()));
            sj.add(topAlignment.getXATag());
            sj.add(String.valueOf(topAlignment.getAlignerScore()));
            sj.add(String.valueOf(nextAlignment != null ? nextAlignment.getAlignerScore() : -1));

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise decoy match writer: {}", e.toString());
        }
    }
}
