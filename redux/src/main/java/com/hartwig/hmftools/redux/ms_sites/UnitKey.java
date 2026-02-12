package com.hartwig.hmftools.redux.ms_sites;

public enum UnitKey
{
    _3_TO_5_BP("3-5pb"),
    A_T("A/T"),
    C_G("C/G"),
    AT_TA("AT/TA"),
    AG_GA_CT_TC("AG/GA/CT/TC"),
    AC_CA_GT_TG("AC/CA/GT/TG"),
    CG_GC("CG/GC");

    private final String unitKey;
    public final int length;

    public String getUnitKey()
    {
        return unitKey;
    }

    UnitKey(String unitKey)
    {
        this.unitKey = unitKey;
        if(unitKey.contains("/"))
        {
            this.length = unitKey.split("/")[0].length();
        }
        else
        {
            this.length = 3;
        }
    }

    public static UnitKey fromUnit(String unit)
    {
        if(unit.length() >= 3 && unit.length() <= 5)
        {
            return _3_TO_5_BP;
        }
        if(unit.equals("A") || unit.equals("T"))
        {
            return A_T;
        }
        if(unit.equals("C") || unit.equals("G"))
        {
            return C_G;
        }
        if(unit.equals("AT") || unit.equals("TA"))
        {
            return AT_TA;
        }
        if(unit.equals("AG") || unit.equals("GA") || unit.equals("CT") || unit.equals("TC"))
        {
            return AG_GA_CT_TC;
        }
        if(unit.equals("AC") || unit.equals("CA") || unit.equals("GT") || unit.equals("TG"))
        {
            return AC_CA_GT_TG;
        }
        if(unit.equals("CG") || unit.equals("GC"))
        {
            return CG_GC;
        }

        throw new IllegalStateException("unable to convert unit key: " + unit);
    }
}
