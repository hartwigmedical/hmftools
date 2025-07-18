package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.READ_ID_TRIMMER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import static htsjdk.samtools.SAMFlag.FIRST_OF_PAIR;
import static htsjdk.samtools.SAMFlag.MATE_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;
import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.READ_UNMAPPED;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class AssemblyReadWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;

    public AssemblyReadWriter(final AssemblyConfig config)
    {
        mConfig = config;
        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.ASSEMBLY_READ))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.ASSEMBLY_READ));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("AssemblyId");
            sj.add("AssemblyInfo"); // for now junction info, in time phased assembly info

            sj.add("ReadId");
            sj.add("SupportType");
            sj.add("IsRef");
            sj.add("Chromosome");
            sj.add("PosStart");
            sj.add("PosEnd");
            sj.add("Orient");
            sj.add("Cigar");
            sj.add("InsertSize");
            sj.add("MateChr");
            sj.add("MatePosStart");
            sj.add("MatePosEnd");
            sj.add("Flags");
            sj.add("FirstInPair");
            sj.add("Reversed");
            sj.add("Unmapped");
            sj.add("MateUnmapped");
            sj.add("MateReversed");
            sj.add("Supplementary");
            sj.add("SuppData");

            sj.add("InferredFragLength");
            sj.add("JunctionReadStartDistance");
            sj.add("AlignmentIndex");
            sj.add("AlignmentOrientation");
            sj.add("BreakendSupport");

            sj.add("Matches");
            sj.add("ExtMismatches");
            sj.add("RefMismatches");
            sj.add("TrimCount");
            sj.add("LineTail");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise assembly writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeAssemblyReads(final JunctionAssembly assembly)
    {
        if(mWriter == null)
            return;

        try
        {
            String assemblyInfo = format("%s", assembly.junction().coords());

            List<SupportRead> supportReads;

            if(AssemblyConfig.WriteCandidateReads)
            {
                supportReads = Lists.newArrayList(assembly.support());
                assembly.candidateSupport().forEach(x -> supportReads.add(
                        new SupportRead(x, SupportType.DISCORDANT, 0, -1, -1)));

                assembly.unmappedReads().forEach(x -> supportReads.add(
                        new SupportRead(x, SupportType.DISCORDANT, 0, -1, -1)));
            }
            else
            {
                supportReads = assembly.support();
            }

            for(SupportRead support : supportReads)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);

                sj.add(String.valueOf(assembly.id()));
                sj.add(assemblyInfo);

                sj.add(READ_ID_TRIMMER.restore(support.id()));
                sj.add(support.type().toString());
                sj.add(String.valueOf(support.isReference()));
                sj.add(support.chromosome());
                sj.add(String.valueOf(support.alignmentStart()));
                sj.add(String.valueOf(support.alignmentEnd()));
                sj.add(String.valueOf(support.orientation().asByte()));
                sj.add(support.cigar());
                sj.add(String.valueOf(support.insertSize()));

                sj.add(support.mateChromosome());
                sj.add(String.valueOf(support.mateAlignmentStart()));
                sj.add(String.valueOf(support.mateAlignmentEnd()));

                sj.add(String.valueOf(support.flags()));
                sj.add(String.valueOf(support.isFlagSet(FIRST_OF_PAIR)));
                sj.add(String.valueOf(support.isFlagSet(READ_REVERSE_STRAND)));
                sj.add(String.valueOf(support.isFlagSet(READ_UNMAPPED)));
                sj.add(String.valueOf(support.isFlagSet(MATE_UNMAPPED)));
                sj.add(String.valueOf(support.isFlagSet(MATE_REVERSE_STRAND)));

                sj.add(String.valueOf(support.isSupplementary()));

                if(support.supplementaryData() != null)
                {
                    sj.add(format("%s:%d:%s",
                            support.supplementaryData().Chromosome, support.supplementaryData().Position, support.supplementaryData().Cigar));
                }
                else
                {
                    sj.add("");
                }

                sj.add(String.valueOf(support.inferredFragmentLength()));
                sj.add(String.valueOf(support.junctionReadStartDistance()));
                sj.add(String.valueOf(support.fullAssemblyIndexStart()));
                sj.add(String.valueOf(support.fullAssemblyOrientation() != null ? support.fullAssemblyOrientation().asByte() : 0));

                sj.add(support.breakendSupportType() != null ? support.breakendSupportType().toString() : "NONE");

                sj.add(String.valueOf(support.extensionBaseMatches()));
                sj.add(String.valueOf(support.extensionBaseMismatches()));
                sj.add(String.valueOf(support.referenceMismatches()));
                sj.add(String.valueOf(support.trimCount()));
                sj.add(String.valueOf(support.hasLineTail()));

                mWriter.write(sj.toString());
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly reads: {}", e.toString());
        }
    }
}
