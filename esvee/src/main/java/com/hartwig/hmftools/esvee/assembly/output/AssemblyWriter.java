package com.hartwig.hmftools.esvee.assembly.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.READ_ID_TRIMMER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.addPhasingHeader;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.addPhasingInfo;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.addRemoteRegionHeader;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.addRemoteRegionInfo;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.addSupportCounts;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.addSupportHeader;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.refSideSoftClipsStr;
import static com.hartwig.hmftools.esvee.assembly.output.AssemblyWriterUtils.repeatsInfoStr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.AssemblyUtils;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyStats;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

public class AssemblyWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;

    // write info about assemblies
    public AssemblyWriter(final AssemblyConfig config)
    {
        mConfig = config;
        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.JUNC_ASSEMBLY))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.JUNC_ASSEMBLY));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Id");
            sj.add("Chromosome").add("JuncPosition").add("JuncOrientation").add("JuncType");

            sj.add("ExtBaseLength").add("RefBasePosition").add("RefBaseLength").add("RefBaseCigar");

            addSupportHeader(sj);

            sj.add("Outcome");

            addPhasingHeader(sj);

            sj.add("JunctionSequence");
            sj.add("RefBaseSequence");
            sj.add("InsertType");

            sj.add("RefBaseCandidates");
            sj.add("UnmappedCandidates");

            sj.add("AssemblyInfo");

            // extra detailed fields
            if(mConfig.AssemblyDetailedTsv)
            {
                sj.add("InitialReadId");
                sj.add("ExtBaseBuildInfo");
                sj.add("MismatchReads");

                sj.add("RefSideSoftClips");
                sj.add("RefBaseTrimmed");
                sj.add("RefBaseTrimLength");
                sj.add("RepeatInfo");
                AssemblyStats.addReadStatsHeader(sj);

                AssemblyStats.addReadTypeHeader(sj);
                addRemoteRegionHeader(sj);
                sj.add("MergedAssemblies");
            }

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

    public void writeAssembly(final JunctionAssembly assembly)
    {
        if(mWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(String.valueOf(assembly.id()));
            sj.add(assembly.junction().Chromosome);
            sj.add(String.valueOf(assembly.junction().Position));
            sj.add(String.valueOf(assembly.junction().Orient));

            String juncType = assembly.junction().DiscordantOnly ? "DISC" :
                    (assembly.junction().indelBased() ? "INDEL" : "SPLIT");
            sj.add(juncType);

            sj.add(String.valueOf(assembly.extensionLength()));
            sj.add(String.valueOf(assembly.refBasePosition()));
            sj.add(String.valueOf(assembly.refBaseLength()));
            sj.add(String.valueOf(assembly.refBaseCigar()));

            addSupportCounts(assembly, sj);

            sj.add(String.valueOf(assembly.outcome()));

            addPhasingInfo(assembly, sj);

            if(AssemblyUtils.hasUnsetBases(assembly))
            {
                sj.add("UNSET_BASES");
                sj.add("UNSET_BASES");
            }
            else
            {
                sj.add(assembly.formJunctionSequence());

                int refBaseLength = mConfig.AssemblyRefBaseWriteMax == 0 ? assembly.refBaseLength() : mConfig.AssemblyRefBaseWriteMax;
                sj.add(assembly.formRefBaseSequence(refBaseLength)); // long enough to show most short TIs
            }

            String insertionType =  assembly.hasLineSequence() ? "LINE" : "NONE";
            sj.add(insertionType);

            sj.add(String.valueOf(assembly.stats().CandidateSupportCount));
            sj.add(String.valueOf(assembly.stats().UnmappedReadCount));

            sj.add(assembly.assemblyAlignmentInfo());

            if(mConfig.AssemblyDetailedTsv)
            {
                sj.add(READ_ID_TRIMMER.restore(assembly.initialReadId()));
                sj.add(assembly.extBaseBuildInfo());
                sj.add(String.valueOf(assembly.mismatchReadCount()));

                sj.add(refSideSoftClipsStr(assembly.refSideSoftClips()));

                sj.add(assembly.refBasesRepeatedTrimmed());
                sj.add(String.valueOf(assembly.refBaseTrimLength()));

                sj.add(repeatsInfoStr(assembly.repeatInfo()));

                assembly.stats().addReadStats(sj);
                assembly.stats().addReadTypeCounts(sj);

                addRemoteRegionInfo(assembly, sj);
                sj.add(String.valueOf(assembly.mergedAssemblyCount()));
            }

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly: {}", e.toString());
        }
    }
}
