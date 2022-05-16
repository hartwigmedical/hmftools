package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class DriverCatalogMap {

    private DriverCatalogMap() {
    }

    @NotNull
    public static Map<DriverCatalogKey, DriverCatalog> toDriverMap(@NotNull List<DriverCatalog> driverCatalog) {
        Map<DriverCatalogKey, DriverCatalog> map = Maps.newHashMap();
        for (DriverCatalog driver : driverCatalog) {
            DriverCatalogKey key = DriverCatalogKey.create(driver.gene(), driver.transcript());
            map.put(key, driver);
        }
        return map;
    }
}
