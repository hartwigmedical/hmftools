package com.hartwig.hmftools.common.utils.config;

import com.beust.jcommander.IStringConverter;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

// we need to define a converter for ref genome version
public class RefGenomeVersionConverter implements IStringConverter<RefGenomeVersion>
{
    @Override
    public RefGenomeVersion convert(String value)
    {
        return RefGenomeVersion.from(value);
    }
}
