package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.DECOY_MAX_MISMATCHES;

import static htsjdk.samtools.CigarOperator.M;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

import htsjdk.samtools.Cigar;

public class AlignmentWriter
{
    private final BwaMemAligner mAligner;
    private final BufferedWriter mWriter;
    private final BufferedWriter mDetailedWriter;

    public AlignmentWriter(final AssemblyConfig config)
    {
        mWriter = initialiseWriter(config);
        mDetailedWriter = initialiseDetailedWriter(config);
        mAligner = null;
    }

    public BufferedWriter alignmentWriter() { return mWriter; }
    public BufferedWriter alignmentDetailsWriter() { return mDetailedWriter; }

    private BufferedWriter initialiseWriter(final AssemblyConfig config)
    {
        if(!config.WriteTypes.contains(WriteType.ALIGNMENT))
            return null;

        if(config.OutputDir == null)
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(config.outputFilename(WriteType.ALIGNMENT));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("AssemblyIds");
            sj.add("AssemblyInfo");
            sj.add("SvType");
            sj.add("SvLength");
            sj.add("SequenceLength");
            sj.add("AssemblyCigar");

            sj.add("AlignCigar");
            sj.add("AlignScore");
            sj.add("AlignedBases");

            sj.add("SecondCigar");
            sj.add("SecondScore");
            sj.add("SecondAlignedBases");

            sj.add("FullSequence");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise alignment writer: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeAssemblyAlignment(
            final BufferedWriter writer, final AssemblyAlignment assemblyAlignment, final String fullSequence,
            final List<BwaMemAlignment> alignmentResults)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(assemblyAlignment.ids());
            sj.add(assemblyAlignment.info());
            sj.add(String.valueOf(assemblyAlignment.svType()));
            sj.add(String.valueOf(assemblyAlignment.svLength()));
            sj.add(String.valueOf(fullSequence.length()));
            sj.add(assemblyAlignment.assemblyCigar());

            BwaMemAlignment topAlignment = !alignmentResults.isEmpty() ? alignmentResults.get(0) : null;

            if(topAlignment == null || topAlignment.getAlignerScore() == 0 || topAlignment.getCigar().isEmpty())
            {
                sj.add("").add("0").add("0").add("").add("0").add("0");
                sj.add(fullSequence);
                writer.write(sj.toString());
                writer.newLine();
                return;
            }

            Cigar cigar = CigarUtils.cigarFromStr(topAlignment.getCigar());
            int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

            sj.add(topAlignment.getCigar());
            sj.add(String.valueOf(topAlignment.getAlignerScore()));
            sj.add(String.valueOf(alignedBases));

            BwaMemAlignment nextAlignment = alignmentResults.size() > 1 ? alignmentResults.get(1) : null;

            if(nextAlignment != null)
            {
                cigar = CigarUtils.cigarFromStr(nextAlignment.getCigar());
                alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

                sj.add(nextAlignment.getCigar());
                sj.add(String.valueOf(nextAlignment.getAlignerScore()));
                sj.add(String.valueOf(alignedBases));
            }
            else
            {
                sj.add("").add("0").add("0");
            }

            sj.add(fullSequence);

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write alignment data: {}", e.toString());
        }
    }

    private BufferedWriter initialiseDetailedWriter(final AssemblyConfig config)
    {
        if(!config.WriteTypes.contains(WriteType.ALIGNMENT_DETAILED))
            return null;

        if(config.OutputDir == null)
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(config.outputFilename(WriteType.ALIGNMENT_DETAILED));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("AssemblyIds");
            sj.add("AssemblyInfo");

            sj.add("RefInfo");
            sj.add("SequenceCoords");
            sj.add("MapQual");
            sj.add("Cigar");
            sj.add("AlignedBases");
            sj.add("Score");
            sj.add("NMatches");
            // sj.add("MateRefInfo");
            // sj.add("TemplateLength");
            sj.add("Location");
            sj.add("MDTag");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise detailed alignment data writer: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeAlignmentDetails(
            final BufferedWriter writer, final AssemblyAlignment assemblyAlignment, final List<BwaMemAlignment> alignments)
    {
        if(writer == null)
            return;

        try
        {
            String assemblyStr = format("%s\t%s", assemblyAlignment.ids(), assemblyAlignment.info());

            for(BwaMemAlignment alignment : alignments)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(assemblyStr);
                sj.add(format("%d:%d-%d", alignment.getRefId(), alignment.getRefStart(), alignment.getRefEnd()));
                sj.add(format("%d-%d", alignment.getSeqStart(), alignment.getSeqEnd()));
                sj.add(String.valueOf(alignment.getMapQual()));
                sj.add(String.valueOf(alignment.getCigar()));

                Cigar cigar = CigarUtils.cigarFromStr(alignment.getCigar());
                int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();
                sj.add(String.valueOf(alignedBases));

                sj.add(String.valueOf(alignment.getAlignerScore()));
                sj.add(String.valueOf(alignment.getNMismatches()));
                // sj.add(format("%d-%d", alignment.getMateRefId(), alignment.getMateRefStart()));
                //sj.add(String.valueOf(alignment.getTemplateLen()));
                sj.add(alignment.getXATag());
                sj.add(alignment.getMDTag());

                writer.write(sj.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write alignment details: {}", e.toString());
        }
    }
}
