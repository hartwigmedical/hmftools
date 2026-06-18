package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

public class ReadJunctionInfo
{
    public final com.hartwig.hmftools.esvee.assembly.read.Read Read;
    public final int ExtensionLength;
    public final boolean MatchesJunction;

    public ReadJunctionInfo(
            final com.hartwig.hmftools.esvee.assembly.read.Read read, final int extensionLength, final boolean matchesJunction)
    {
        Read = read;
        ExtensionLength = extensionLength;
        MatchesJunction = matchesJunction;
    }

    public String toString()
    {
        return format("%s: extension(%s) juncMatch(%s)", Read.id(), ExtensionLength, MatchesJunction);
    }
}
