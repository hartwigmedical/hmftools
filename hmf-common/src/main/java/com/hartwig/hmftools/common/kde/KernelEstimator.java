package com.hartwig.hmftools.common.kde;

// JOBA: This is a modified version of the WEKA KernelEstimator
public class KernelEstimator   {
    private double[] m_Values = new double[50];
    private double[] m_Weights = new double[50];
    private int m_NumValues = 0;
    private double m_SumOfWeights = 0.0D;
    private final double m_StandardDev;
    private double m_Precision;
    private static final double MAX_ERROR = 0.01D;

    public KernelEstimator(double precision, double bandwidth) {
        this.m_Precision = precision;
        this.m_StandardDev = bandwidth;
    }

    public void addValue(double data, double weight) {
        if(weight != 0.0D) {
            data = this.round(data);
            int insertIndex = this.findNearestValue(data);
            if(this.m_NumValues > insertIndex && this.m_Values[insertIndex] == data) {
                this.m_Weights[insertIndex] += weight;
            } else {
                if(this.m_NumValues < this.m_Values.length) {
                    int left = this.m_NumValues - insertIndex;
                    System.arraycopy(this.m_Values, insertIndex, this.m_Values, insertIndex + 1, left);
                    System.arraycopy(this.m_Weights, insertIndex, this.m_Weights, insertIndex + 1, left);
                    this.m_Values[insertIndex] = data;
                    this.m_Weights[insertIndex] = weight;
                    ++this.m_NumValues;
                } else {
                    double[] newValues = new double[this.m_Values.length * 2];
                    double[] newWeights = new double[this.m_Values.length * 2];
                    int left = this.m_NumValues - insertIndex;
                    System.arraycopy(this.m_Values, 0, newValues, 0, insertIndex);
                    System.arraycopy(this.m_Weights, 0, newWeights, 0, insertIndex);
                    newValues[insertIndex] = data;
                    newWeights[insertIndex] = weight;
                    System.arraycopy(this.m_Values, insertIndex, newValues, insertIndex + 1, left);
                    System.arraycopy(this.m_Weights, insertIndex, newWeights, insertIndex + 1, left);
                    ++this.m_NumValues;
                    this.m_Values = newValues;
                    this.m_Weights = newWeights;
                }

            }

            this.m_SumOfWeights += weight;
        }
    }

    public double getProbability(double data) {
        double delta;
        double sum = 0.0D;
        double currentProb;
        double zLower;
        double zUpper;
        if(this.m_NumValues == 0) {
            zLower = (data - this.m_Precision / 2.0D) / this.m_StandardDev;
            zUpper = (data + this.m_Precision / 2.0D) / this.m_StandardDev;
            return Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
        } else {
            double weightSum = 0.0D;
            int start = this.findNearestValue(data);

            int i;
            for(i = start; i < this.m_NumValues; ++i) {
                delta = this.m_Values[i] - data;
                zLower = (delta - this.m_Precision / 2.0D) / this.m_StandardDev;
                zUpper = (delta + this.m_Precision / 2.0D) / this.m_StandardDev;
                currentProb = Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
                sum += currentProb * this.m_Weights[i];
                weightSum += this.m_Weights[i];
                if(currentProb * (this.m_SumOfWeights - weightSum) < sum * MAX_ERROR) {
                    break;
                }
            }

            for(i = start - 1; i >= 0; --i) {
                delta = this.m_Values[i] - data;
                zLower = (delta - this.m_Precision / 2.0D) / this.m_StandardDev;
                zUpper = (delta + this.m_Precision / 2.0D) / this.m_StandardDev;
                currentProb = Statistics.normalProbability(zUpper) - Statistics.normalProbability(zLower);
                sum += currentProb * this.m_Weights[i];
                weightSum += this.m_Weights[i];
                if(currentProb * (this.m_SumOfWeights - weightSum) < sum * MAX_ERROR) {
                    break;
                }
            }

            return sum / this.m_SumOfWeights;
        }
    }

    private int findNearestValue(double key) {
        int low = 0;
        int high = this.m_NumValues;

        while(low < high) {
            int middle = (low + high) / 2;
            double current = this.m_Values[middle];
            if(current == key) {
                return middle;
            }

            if(current > key) {
                high = middle;
            } else if(current < key) {
                low = middle + 1;
            }
        }

        return low;
    }

    private double round(double data) {
        return Math.rint(data / this.m_Precision) * this.m_Precision;
    }
}
