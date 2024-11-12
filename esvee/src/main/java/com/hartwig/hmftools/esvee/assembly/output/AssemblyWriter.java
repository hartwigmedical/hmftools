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
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

public class AssemblyWriter
{
    private final AssemblyConfig mConfig;

    private final BufferedWriter mWriter;
    private final TruthsetAnnotation mTruthsetAnnotation;

    // write info about assemblies
    public AssemblyWriter(final AssemblyConfig config, final TruthsetAnnotation truthsetAnnotation)
    {
        mConfig = config;
        mTruthsetAnnotation = truthsetAnnotation;

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

            sj.add("MismatchReads");

            sj.add("RefBaseTrimmed");
            sj.add("RefBaseTrimLength");
            sj.add("JunctionSequence");
            sj.add("RefBaseSequence");
            sj.add("IsLINE");

            sj.add("RefBaseCandidates");
            sj.add("UnmappedCandidates");
            AssemblyStats.addReadTypeHeader(sj);
            addRemoteRegionHeader(sj);

            sj.add("AssemblyInfo");

            // extra detailed fields
            sj.add("InitialReadId");

            sj.add("RepeatInfo");
            sj.add("RefSideSoftClips");
            sj.add("MergedAssemblies");
            AssemblyStats.addReadStatsHeader(sj);
            sj.add("ExtBaseBuildInfo");

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

            sj.add(String.valueOf(assembly.mismatchReadCount()));

            sj.add(assembly.refBasesRepeatedTrimmed());
            sj.add(String.valueOf(assembly.refBaseTrimLength()));

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

            sj.add(String.valueOf(assembly.hasLineSequence()));
            sj.add(String.valueOf(assembly.stats().CandidateSupportCount));
            sj.add(String.valueOf(assembly.stats().UnmappedReadCount));
            assembly.stats().addReadTypeCounts(sj);

            addRemoteRegionInfo(assembly, sj);

            sj.add(assembly.assemblyAlignmentInfo());

            sj.add(READ_ID_TRIMMER.restore(assembly.initialReadId()));

            sj.add(repeatsInfoStr(assembly.repeatInfo()));

            sj.add(refSideSoftClipsStr(assembly.refSideSoftClips()));
            sj.add(String.valueOf(assembly.mergedAssemblyCount()));
            assembly.stats().addReadStats(sj);
            sj.add(assembly.extBaseBuildInfo());

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly: {}", e.toString());
        }
    }
}
