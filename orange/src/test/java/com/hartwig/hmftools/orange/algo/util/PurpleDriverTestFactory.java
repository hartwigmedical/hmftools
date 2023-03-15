package com.hartwig.hmftools.orange.algo.util;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpreter;

public class PurpleDriverTestFactory {
    public static ImmutablePurpleDriver.Builder builder() {
        var builder = DriverCatalogTestFactory.builder();
        var driver = PurpleInterpreter.asPurpleDriver(builder.build());
        return ImmutablePurpleDriver.builder().from(driver);
    }
}
