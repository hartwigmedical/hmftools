package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;

import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.types.AssemblyLink;
import com.hartwig.hmftools.esvee.types.PhaseSet;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class Alignment
{
    private final AssemblyConfig mConfig;

    private final Aligner mAligner;

    public Alignment(final AssemblyConfig config, final Aligner aligner)
    {
        mConfig = config;
        mAligner = aligner;
    }

    public void processPhaseSet(final PhaseSet phaseSet)
    {
        final byte[] fullAssembly = phaseSet.buildFullAssembly();

        if(fullAssembly == null)
            return;

        AssemblyLink assemblyLink = phaseSet.assemblyLinks().get(0);
        List<BwaMemAlignment> alignments = mAligner.alignSequence(fullAssembly);

        for(BwaMemAlignment alignment : alignments)
        {
            SV_LOGGER.debug("assemblyLink({}) alignment: {}", assemblyLink, alignment);
        }

    }

    public static String alignmentStr(final BwaMemAlignment alignment)
    {
        return format("");
    }


}
