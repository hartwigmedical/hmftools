package com.hartwig.hmftools.common.drivercatalog;

import java.util.EnumSet;

public enum DriverType
{
    AMP,
    PARTIAL_AMP,
    DEL,
    HOM_DUP_DISRUPTION,
    HOM_DEL_DISRUPTION,
    MUTATION,
    GERMLINE, // to be deprecated, was previously just germline variants
    GERMLINE_MUTATION,
    GERMLINE_DELETION,
    GERMLINE_DISRUPTION;

    public static final EnumSet<DriverType> DRIVERS_PURPLE_GERMLINE = EnumSet.of(DriverType.GERMLINE_MUTATION, DriverType.GERMLINE_DELETION);

    public static final EnumSet<DriverType> DRIVERS_PURPLE_SOMATIC = EnumSet.of(
            DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL, DriverType.MUTATION);

    public static final EnumSet<DriverType> DRIVERS_LINX_GERMLINE = EnumSet.of(DriverType.GERMLINE_DISRUPTION);

    public static final EnumSet<DriverType> DRIVERS_LINX_SOMATIC = EnumSet.of(DriverType.HOM_DUP_DISRUPTION, DriverType.HOM_DEL_DISRUPTION);

    public static boolean isGermline(final DriverType type)
    {
        return type == GERMLINE_DELETION || type == GERMLINE_DISRUPTION || type == GERMLINE_MUTATION || type == GERMLINE;
    }

    public static DriverType checkConvertType(final String driverTypeStr)
    {
        // support old Purple and Linx names prior to splits for new germline types
        if(driverTypeStr.equals("HOM_DISRUPTION"))
            return HOM_DEL_DISRUPTION;

        if(driverTypeStr.equals("GERMLINE"))
            return GERMLINE_MUTATION;

        return DriverType.valueOf(driverTypeStr);
    }

}
