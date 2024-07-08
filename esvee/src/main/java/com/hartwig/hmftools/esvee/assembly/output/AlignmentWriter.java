package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import static htsjdk.samtools.CigarOperator.M;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;

import htsjdk.samtools.Cigar;

public class AlignmentWriter
{
    private final BufferedWriter mWriter;
    private final BufferedWriter mDetailedWriter;

    public AlignmentWriter(final AssemblyConfig config)
    {
        if(config.AlignmentFile == null)
        {
            mWriter = initialiseWriter(config);
            mDetailedWriter = initialiseDetailedWriter(config);
        }
        else
        {
            mWriter = null;
            mDetailedWriter = null;
        }
    }

    public BufferedWriter alignmentWriter() { return mWriter; }
    public BufferedWriter alignmentDetailsWriter() { return mDetailedWriter; }

    public void close()
    {
        closeBufferedWriter(mWriter);
        closeBufferedWriter(mDetailedWriter);
    }

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

            sj.add(FLD_ASSEMBLY_IDS);
            sj.add(FLD_ASSEMLY_INFO);
            sj.add("SvType");
            sj.add("SvLength");
            sj.add("RefBaseLength");
            sj.add("SequenceLength");
            sj.add("AssemblyCigar");

            sj.add("AlignResults");
            sj.add("AlignCigar");
            sj.add("AlignScore");
            sj.add("AlignedBases");

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
            final BufferedWriter writer, final AssemblyAlignment assemblyAlignment, final List<AlignData> alignmentResults)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(assemblyAlignment.assemblyIds());
            sj.add(assemblyAlignment.info());
            sj.add(String.valueOf(assemblyAlignment.svType()));
            sj.add(String.valueOf(assemblyAlignment.svLength()));
            sj.add(String.valueOf(assemblyAlignment.refBaseLength()));
            sj.add(String.valueOf(assemblyAlignment.fullSequenceLength()));
            sj.add(assemblyAlignment.assemblyCigar());

            AlignData topAlignment = !alignmentResults.isEmpty() ? alignmentResults.get(0) : null;

            if(topAlignment == null || topAlignment.Score == 0 || topAlignment.Cigar.isEmpty())
            {
                sj.add("0").add("").add("0").add("0").add("").add("0");
                sj.add(assemblyAlignment.fullSequence());
                writer.write(sj.toString());
                writer.newLine();
                return;
            }

            sj.add(String.valueOf(alignmentResults.size()));

            Cigar cigar = CigarUtils.cigarFromStr(topAlignment.Cigar);
            int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

            sj.add(topAlignment.Cigar);
            sj.add(String.valueOf(topAlignment.Score));
            sj.add(String.valueOf(alignedBases));

            sj.add(assemblyAlignment.fullSequence());

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write alignment data: {}", e.toString());
        }
    }

    public static final String FLD_ASSEMBLY_IDS = "AssemblyIds";
    public static final String FLD_ASSEMLY_INFO = "AssemblyInfo";
    public static final String FLD_REF_LOCATION = "RefInfo";
    public static final String FLD_SEQUENCE_COORDS = "SequenceCoords";
    public static final String FLD_MAP_QUAL = "MapQual";
    public static final String FLD_CIGAR = "Cigar";
    public static final String FLD_ALIGNED_BASES = "AlignedBases";
    public static final String FLD_SCORE = "Score";
    public static final String FLD_FLAGS = "Flags";
    public static final String FLD_NMATCHES = "NMatches";
    public static final String FLD_XA_TAG = "LocTag";
    public static final String FLD_MD_TAG = "MdTag";

    private BufferedWriter initialiseDetailedWriter(final AssemblyConfig config)
    {
        if(!config.WriteTypes.contains(WriteType.ALIGNMENT_DATA))
            return null;

        if(config.OutputDir == null)
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(config.outputFilename(WriteType.ALIGNMENT_DATA));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_ASSEMBLY_IDS);
            sj.add(FLD_ASSEMLY_INFO);

            sj.add(FLD_REF_LOCATION);
            sj.add(FLD_SEQUENCE_COORDS);
            sj.add(FLD_MAP_QUAL);
            sj.add(FLD_CIGAR);
            sj.add(FLD_ORIENTATION);
            sj.add(FLD_ALIGNED_BASES);
            sj.add(FLD_SCORE);
            sj.add(FLD_FLAGS);
            sj.add(FLD_NMATCHES);
            sj.add(FLD_XA_TAG);
            sj.add(FLD_MD_TAG);
            sj.add("Requeried");

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
            final BufferedWriter writer, final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments)
    {
        if(writer == null || alignments.isEmpty())
            return;

        try
        {
            String assemblyStr = format("%s\t%s", assemblyAlignment.assemblyIds(), assemblyAlignment.info());

            for(AlignData alignment : alignments)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(assemblyStr);
                sj.add(alignment.RefLocation.toString());
                sj.add(format("%d-%d", alignment.rawSequenceStart(), alignment.rawSequenceEnd()));
                sj.add(String.valueOf(alignment.MapQual));
                sj.add(String.valueOf(alignment.Cigar));
                sj.add(String.valueOf(alignment.orientation()));

                Cigar cigar = CigarUtils.cigarFromStr(alignment.Cigar);
                int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();
                sj.add(String.valueOf(alignedBases));

                sj.add(String.valueOf(alignment.Score));
                sj.add(String.valueOf(alignment.Flags));
                sj.add(String.valueOf(alignment.NMatches));
                sj.add(alignment.XaTag);
                sj.add(alignment.MdTag);
                sj.add(String.valueOf(alignment.isRequeried()));

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
