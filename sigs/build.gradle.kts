
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Sigs"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.sigs.fitter.SampleFitter")
}
