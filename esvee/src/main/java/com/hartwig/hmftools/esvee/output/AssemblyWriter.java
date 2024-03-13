package com.hartwig.hmftools.esvee.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addPhasingHeader;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addPhasingInfo;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addReadStats;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addReadStatsHeader;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addRemoteRegionHeader;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addRemoteRegionInfo;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addSupportCounts;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.addSupportHeader;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.refSideSoftClipsStr;
import static com.hartwig.hmftools.esvee.output.AssemblyWriterUtils.repeatsInfoStr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

public class AssemblyWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mWriter;
    private final TruthsetAnnotation mTruthsetAnnotation;

    // write info about assemblies
    public AssemblyWriter(final SvConfig config)
    {
        mConfig = config;
        mTruthsetAnnotation = new TruthsetAnnotation(mConfig.TruthsetFile);

        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter);}

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.ASSEMBLIES))
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.ASSEMBLIES));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Id");
            sj.add("Chromosome").add("JunctionPosition").add("JunctionOrientation");

            sj.add("ExtBaseLength").add("RefBasePosition").add("RefBaseLength");

            addSupportHeader(sj);

            addReadStatsHeader(sj);

            sj.add("Outcome");

            addPhasingHeader(sj);

            if(mTruthsetAnnotation.enabled())
                sj.add(TruthsetAnnotation.tsvHeader());

            addRemoteRegionHeader(sj);

            sj.add("RefBaseTrimmed");
            sj.add("JunctionSequence");
            sj.add("RefBaseSequence");

            // extra detailed fields
            sj.add("InitialReadId");

            sj.add("InitRefBaseCandidates");

            sj.add("MergedAssemblies");

            sj.add("RepeatInfo");
            sj.add("RefSideSoftClips");
            sj.add("BranchedAssemblyIds");

            if(mConfig.LogPhaseGroupLinks)
                sj.add("PhaseGroupLinkInfo");

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
            sj.add(String.valueOf(assembly.junction().Orientation));

            sj.add(String.valueOf(assembly.extensionLength()));
            sj.add(String.valueOf(assembly.refBasePosition()));
            sj.add(String.valueOf(assembly.refBaseLength()));

            addSupportCounts(assembly, sj);

            addReadStats(assembly, sj);

            sj.add(String.valueOf(assembly.outcome()));

            addPhasingInfo(assembly, sj);

            if(mTruthsetAnnotation.enabled())
                sj.add(mTruthsetAnnotation.findTruthsetAnnotation(assembly));

            addRemoteRegionInfo(assembly, sj);

            sj.add(assembly.refBasesRepeatedTrimmed());

            if(assembly.hasUnsetBases())
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

            sj.add(assembly.initialRead() != null ? assembly.initialRead().getName() : "NONE"); // shouldn't occur

            sj.add(String.valueOf(assembly.candidateSupport().size()));

            sj.add(String.valueOf(assembly.mergedAssemblyCount()));

            sj.add(repeatsInfoStr(assembly.repeatInfo()));

            sj.add(refSideSoftClipsStr(assembly.refSideSoftClips()));

            String branchedAssemblyIds = assembly.branchedAssemblies().stream()
                    .map(x -> String.valueOf(x.id())).collect(Collectors.joining(";"));
            sj.add(branchedAssemblyIds);

            if(mConfig.LogPhaseGroupLinks)
                sj.add(assembly.phaseGroupLinkingInfo());

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write assembly: {}", e.toString());
        }
    }
}
