
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Isofox"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.isofox.Isofox")
}
