package com.hartwig.hmftools.esvee.assembly.output;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;

import java.io.BufferedWriter;

import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignmentChecker;

public class ResultsWriter
{
    private final BufferedWriter mDecoyMatchWriter;
    private final AssemblyWriter mAssemblyWriter;
    private final BreakendWriter mBreakendWriter;
    private final AssemblyReadWriter mReadWriter;
    private final PhaseGroupBuildWriter mPhaseGroupBuildWriter;
    private final BamWriter mBamWriter;

    public ResultsWriter(final AssemblyConfig config)
    {
        mAssemblyWriter = new AssemblyWriter(config);
        mBreakendWriter = new BreakendWriter(config);
        mReadWriter = new AssemblyReadWriter(config);
        mPhaseGroupBuildWriter = new PhaseGroupBuildWriter(config);
        mBamWriter = new BamWriter(config);
        mDecoyMatchWriter = AlignmentChecker.initialiseWriter(config);
    }

    public BufferedWriter decoyMatchWriter() { return mDecoyMatchWriter; }
    public AssemblyWriter assemblyWriter() { return mAssemblyWriter; }
    public BreakendWriter breakendWriter() { return mBreakendWriter; }
    public AssemblyReadWriter readWriter() { return mReadWriter; }
    public PhaseGroupBuildWriter phaseGroupBuildWriter() { return mPhaseGroupBuildWriter; }
    public BamWriter bamWriter() { return mBamWriter; }

    public void close()
    {
        mAssemblyWriter.close();
        mBreakendWriter.close();
        mReadWriter.close();
        mPhaseGroupBuildWriter.close();
        mBamWriter.close();
        closeBufferedWriter(mDecoyMatchWriter);
    }
}
