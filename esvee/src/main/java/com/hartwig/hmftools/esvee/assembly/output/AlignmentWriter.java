package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;

public class AlignmentWriter
{
    private final BufferedWriter mWriter;
    private final BufferedWriter mDetailedWriter;

    public AlignmentWriter(final AssemblyConfig config)
    {
        if(config.AlignmentFile == null)
        {
            mWriter = initialiseWriter(config);
            mDetailedWriter = initialiseAlignmentDataWriter(config);
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
        if(!config.WriteTypes.contains(WriteType.PHASED_ASSEMBLY))
            return null;

        if(config.OutputDir == null)
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(config.outputFilename(WriteType.PHASED_ASSEMBLY));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_ASSEMBLY_IDS);
            sj.add(FLD_ASSEMLY_INFO);
            sj.add("Merged");
            sj.add("SequenceLength");
            sj.add("AssemblyCigar");
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
            final BufferedWriter writer, final AssemblyAlignment assemblyAlignment)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(assemblyAlignment.assemblyIds());
            sj.add(assemblyAlignment.info());
            sj.add(String.valueOf(assemblyAlignment.isMerged()));
            sj.add(String.valueOf(assemblyAlignment.fullSequenceLength()));
            sj.add(assemblyAlignment.assemblyCigar());
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
    public static final String FLD_RAW_SEQ_COORDS = "RawSeqCoords";
    public static final String FLD_ADJ_SEQ_COORDS = "AdjSeqCoords";
    public static final String FLD_MAP_QUAL = "MapQual";
    public static final String FLD_CIGAR = "Cigar";
    public static final String FLD_ALIGNED_BASES = "AlignedBases";
    public static final String FLD_SCORE = "Score";
    public static final String FLD_FLAGS = "Flags";
    public static final String FLD_NMATCHES = "NMatches";
    public static final String FLD_XA_TAG = "LocTag";
    public static final String FLD_MD_TAG = "MdTag";

    private BufferedWriter initialiseAlignmentDataWriter(final AssemblyConfig config)
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

            sj.add(FLD_REF_LOCATION);
            sj.add(FLD_RAW_SEQ_COORDS);
            sj.add(FLD_ADJ_SEQ_COORDS);
            sj.add(FLD_MAP_QUAL);
            sj.add(FLD_CIGAR);
            sj.add(FLD_ORIENTATION);
            sj.add(FLD_ALIGNED_BASES);
            sj.add(FLD_SCORE);
            sj.add(FLD_FLAGS);
            sj.add(FLD_NMATCHES);
            sj.add(FLD_XA_TAG);
            sj.add(FLD_MD_TAG);
            sj.add("CalcAlignLength");
            sj.add("ModMapQual");
            sj.add("DroppedOnRequery");
            sj.add("LinkedAltAlignment");

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
                sj.add(alignment.refLocation().toString());
                sj.add(format("%d-%d", alignment.rawSequenceStart(), alignment.rawSequenceEnd()));
                sj.add(format("%d-%d", alignment.sequenceStart(), alignment.sequenceEnd()));
                sj.add(String.valueOf(alignment.mapQual()));
                sj.add(String.valueOf(alignment.cigar()));
                sj.add(String.valueOf(alignment.orientation()));

                sj.add(String.valueOf(alignment.alignedBases()));

                sj.add(String.valueOf(alignment.score()));
                sj.add(String.valueOf(alignment.flags()));
                sj.add(String.valueOf(alignment.nMatches()));
                sj.add(alignment.xaTag());
                sj.add(alignment.mdTag());
                sj.add(String.valueOf(alignment.adjustedAlignment()));
                sj.add(format("%.0f", alignment.modifiedMapQual()));
                sj.add(String.valueOf(alignment.droppedOnRequery()));

                if(alignment.hasSelectedAltAlignment())
                    sj.add(alignment.selectedAltAlignment().vcfString());
                else
                    sj.add("");

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
