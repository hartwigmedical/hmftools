package com.hartwig.hmftools.esvee.assembly.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
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

import com.hartwig.hmftools.esvee.AssemblyConfig;
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
            sj.add("Chromosome").add("JunctionPosition").add("JunctionOrientation");

            sj.add("ExtBaseLength").add("RefBasePosition").add("RefBaseLength");

            addSupportHeader(sj);
            AssemblyStats.addReadTypeHeader(sj);

            sj.add("Outcome");

            addPhasingHeader(sj);

            AssemblyStats.addReadStatsHeader(sj);
            sj.add("MismatchReads");

            sj.add("RefBaseTrimmed");
            sj.add("RefBaseTrimLength");
            sj.add("JunctionSequence");
            sj.add("RefBaseSequence");

            addRemoteRegionHeader(sj);

            if(mConfig.RunAlignment)
            {
                sj.add("AlignResult");
                sj.add("AssemblyInfo");
            }

            // extra detailed fields
            sj.add("InitialReadId");

            sj.add("InitRefBaseCandidates");

            sj.add("MergedAssemblies");

            sj.add("RepeatInfo");
            sj.add("RefSideSoftClips");

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

            sj.add(String.valueOf(assembly.extensionLength()));
            sj.add(String.valueOf(assembly.refBasePosition()));
            sj.add(String.valueOf(assembly.refBaseLength()));

            addSupportCounts(assembly, sj);
            assembly.stats().addReadTypeCounts(sj);

            sj.add(String.valueOf(assembly.outcome()));

            addPhasingInfo(assembly, sj);

            assembly.stats().addReadStats(sj);
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

            addRemoteRegionInfo(assembly, sj);

            if(mConfig.RunAlignment)
            {
                sj.add(String.valueOf(assembly.alignmentOutcome()));
                sj.add(assembly.assemblyAlignmentInfo());
            }

            sj.add(assembly.initialReadId());

            sj.add(String.valueOf(assembly.stats().CandidateSupportCount));

            sj.add(String.valueOf(assembly.mergedAssemblyCount()));

            sj.add(repeatsInfoStr(assembly.repeatInfo()));

            sj.add(refSideSoftClipsStr(assembly.refSideSoftClips()));

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly: {}", e.toString());
        }
    }
}
