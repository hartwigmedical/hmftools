package com.hartwig.hmftools.qsee.status;

public enum ComparisonOperator
{
    LESS_THAN("<"),
    LESS_THAN_OR_EQUAL("<="),
    GREATER_THAN_OR_EQUAL(">="),
    GREATER_THAN(">");

    public final String mOperatorString;

    ComparisonOperator(String operatorString)
    {
        mOperatorString = operatorString;
    }

    public String operatorString() { return mOperatorString; }

    public static ComparisonOperator fromString(String string)
    {
        if(string.equals(LESS_THAN.toString()) || string.equals(LESS_THAN.operatorString()))
            return LESS_THAN;

        if(string.equals(LESS_THAN_OR_EQUAL.toString()) || string.equals(LESS_THAN_OR_EQUAL.operatorString()))
            return LESS_THAN_OR_EQUAL;

        if(string.equals(GREATER_THAN_OR_EQUAL.toString()) || string.equals(GREATER_THAN_OR_EQUAL.operatorString()))
            return GREATER_THAN_OR_EQUAL;

        if(string.equals(GREATER_THAN.toString()) || string.equals(GREATER_THAN.operatorString()))
            return GREATER_THAN;

        throw new IllegalArgumentException("Unknown comparison operator: " + string);
    }
}
