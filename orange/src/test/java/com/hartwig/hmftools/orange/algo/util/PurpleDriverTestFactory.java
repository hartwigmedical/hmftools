package com.hartwig.hmftools.orange.algo.util;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.NotNull;

public final class PurpleDriverTestFactory
{
    @NotNull
    public static ImmutablePurpleDriver.Builder builder()
    {
        DriverCatalog catalog = DriverCatalogTestFactory.builder().build();
        PurpleDriver driver = PurpleConversion.convert(catalog);
        return ImmutablePurpleDriver.builder().from(driver);
    }
}
