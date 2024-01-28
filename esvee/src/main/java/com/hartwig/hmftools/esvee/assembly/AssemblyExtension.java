package com.hartwig.hmftools.esvee.assembly;

import java.util.List;

import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.read.Read;

public final class AssemblyExtension
{
    public static void extendAssembly(final JunctionAssembly assembly, final List<Read> nonJunctionReads)
    {
        // first establish potential boundaries for extending the assembly on the non-junction side


        /*
        for(Read read : otherReads)
        {
            for(JunctionAssembly assembly : initialAssemblies)
            {
                if(assembly.checkReadMatches(read, PRIMARY_ASSEMBLY_READ_MAX_BASE_MISMATCH))
                {
                    assembly.addRead(read, false);
                }
            }
        }
        */

    }

}
