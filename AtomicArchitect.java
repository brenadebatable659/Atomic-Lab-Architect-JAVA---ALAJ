/**  

MIT License 
Copyright (c) [2026] Sam78887
Created by: Sam78887
Date: March 2026
(See LICENSE file for details)

*/
/**

Atomic Lab Architect JAVA - ALAJ
 * Atomic Architect
 * A Java-based chemical analysis tool and molecular dynamics simulator.
 * * Features
 * * Electron Configuration
 * Automatically generates subshell notation (e.g., 1s^2 2s^2) for elements 1 through 118.
 * * VSEPR Theory
 * Predicts molecular geometry, bond angles, and steric numbers to determine the 3D shape.
 * * MD Simulation
 * Runs a digital lab simulation using a Lennard-Jones potential to model particle interactions.
 * * Physical State Prediction
 * Estimates boiling points and predicts whether a compound will be a solid, liquid, or gas.
 */
import java.util.*;

public class AtomicArchitect {

    static class Atom {
        String name, symbol, config;
        int valenceElectrons, valency, atomicNum;
        double weight, en, radius;

        static final String[] NAMES = {"", "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Assenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
        static final String[] SYMBOLS = {"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
        static final double[] WEIGHTS = new double[119];
        static final double[] EN_SCALE = new double[119];
        static final double[] RADII = new double[119];

        static {
            Arrays.fill(WEIGHTS, 1.0); Arrays.fill(EN_SCALE, 2.0); Arrays.fill(RADII, 100.0);
            double[] w = {0, 1.00784, 4.0026, 6.941, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.98, 28.08, 30.97, 32.06, 35.45, 39.95, 39.10, 40.08};
            System.arraycopy(w, 0, WEIGHTS, 0, w.length);
            WEIGHTS[53] = 126.90; RADII[1]=37; RADII[6]=77; RADII[8]=73; RADII[53]=133;
            EN_SCALE[1]=2.20; EN_SCALE[6]=2.55; EN_SCALE[8]=3.44; EN_SCALE[53]=2.66;
        }

        Atom(int num) {
            this.atomicNum = (num >= 1 && num <= 118) ? num : 0;
            this.name = NAMES[this.atomicNum];
            this.symbol = SYMBOLS[this.atomicNum];
            this.weight = WEIGHTS[this.atomicNum];
            this.en = EN_SCALE[this.atomicNum];
            this.radius = RADII[this.atomicNum];
            generateConfigAndValence(this.atomicNum);
        }

        private void generateConfigAndValence(int n) {
            String[] subshells = {"1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p"};
            int[] caps = {2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6};
            StringBuilder sb = new StringBuilder();
            int temp = n;
            for (int i = 0; i < subshells.length && temp > 0; i++) {
                int fill = Math.min(temp, caps[i]);
                sb.append(subshells[i]).append("^").append(fill).append(" ");
                temp -= fill;
            }
            this.config = sb.toString().trim();
            if (n==2 || n==10 || n==18 || n==36 || n==54 || n==86 || n==118) {
                this.valency = 0; this.valenceElectrons = 8;
            } else {
                int mod = n % 18;
                this.valenceElectrons = (mod == 0) ? 8 : (mod > 10 ? mod - 10 : (mod <= 2 ? mod : (mod <= 12 ? 2 : mod - 10)));
                if (n == 1) valenceElectrons = 1; 
                if (n == 6) valenceElectrons = 4;
                if (n == 7) valenceElectrons = 5;
                if (n == 8) valenceElectrons = 6;
                if (n == 9) valenceElectrons = 7;
                this.valency = (valenceElectrons <= 4) ? valenceElectrons : 8 - valenceElectrons;
            }
        }
    }

    static class Particle {
        double x, y, z, vx, vy, vz, fx, fy, fz;
        Particle(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        String choice = "y";

        while (choice.equalsIgnoreCase("y")) {
            System.out.print("\nHow many unique element types? ");
            int typeCount = sc.nextInt();
            Map<String, Integer> counts = new LinkedHashMap<>();
            List<Atom> atomList = new ArrayList<>();
            double totalWeight = 0;
            int totalValencySum = 0;
            int maxVal = -1;
            Atom central = null, secondary = null;
            int totalAtoms = 0;

            for (int i = 0; i < typeCount; i++) {
                System.out.print("Enter Atomic Number (1-118): ");
                int z = sc.nextInt();
                System.out.print("Enter Quantity: ");
                int qty = sc.nextInt();
                Atom cur = new Atom(z);
                System.out.println("-> " + cur.name + " [" + cur.config + "]");
                totalWeight += (cur.weight * qty);
                totalValencySum += (cur.valency * qty);
                totalAtoms += qty;
                counts.put(cur.symbol, qty);
                for(int j=0; j<qty; j++) atomList.add(cur);

                if (cur.valency > maxVal || (cur.valency == maxVal && qty == 1)) {
                    maxVal = cur.valency;
                    secondary = central;
                    central = cur;
                } else if (secondary == null && !cur.symbol.equals(central.symbol)) {
                    secondary = cur;
                }
            }

            System.out.println("\n--- Discovery Scanner Analysis ---");
            boolean valCheck = (totalValencySum % 2 == 0);
            boolean saturated = true;
            if (counts.containsKey("C")) {
                int others = totalAtoms - counts.get("C");
                if (others < (counts.get("C") / 2 + 1) && !counts.containsKey("O")) saturated = false;
            }

            if (!valCheck || totalValencySum < 2 * maxVal) System.out.println("Status: [IMPOSSIBLE]");
            else if (!saturated) System.out.println("Status: [HIGHLY REACTIVE / THEORETICAL]");
            else {
                if (secondary != null && secondary.radius > central.radius * 1.5 && counts.get(secondary.symbol) >= 4)
                    System.out.println("Status: [PHYSICALLY STRAINED]");
                else System.out.println("Status: [STABLE COMPOUND]");
            }

            System.out.println("\n--- Mass Composition ---");
            for (String s : counts.keySet()) {
                double atomW = 0;
                for (int i=1; i<119; i++) if (Atom.SYMBOLS[i].equals(s)) atomW = Atom.WEIGHTS[i];
                System.out.printf("%s: %.2f%%\n", s, (atomW * counts.get(s) / totalWeight) * 100);
            }

            predictPhysicalState(counts, totalWeight);

            System.out.println("\n--- Structural Details ---");
            System.out.println("Molecular Formula: " + generateFormula(counts));
            System.out.println("SMILES Notation: " + generateSMILES(central, counts));
            System.out.printf("Total Molecular Weight: %.5f g/mol\n", totalWeight);
            
            if (secondary != null) {
                double diff = Math.abs(central.en - secondary.en);
                String bType = (diff >= 1.7) ? "Ionic" : (diff >= 0.5) ? "Polar Covalent" : "Non-Polar Covalent";
                System.out.println("Primary Bond Type: " + bType + " (Delta EN: " + String.format("%.2f", diff) + ")");
                
                // New: Spectroscopy & Bond Analysis
                double bondLength = (central.radius + secondary.radius) - 9 * diff;
                System.out.printf("Predicted Bond Length (%s-%s): %.2f pm\n", central.symbol, secondary.symbol, bondLength);
                double mu = (central.weight * secondary.weight) / (central.weight + secondary.weight);
                //double vibFreq = (1.0 / (2 * Math.PI)) * Math.sqrt((500.0 * (central.en * secondary.en / 4.0)) / (mu * 1.66e-27)) * 1e-12 / 1e8;
                //System.out.printf("Vibration Frequency: ~%.2f THz\n", vibFreq);
                if (secondary != null) {
    // 1. Calculate Reduced Mass (mu) in amu
    double muAmu = (central.weight * secondary.weight) / (central.weight + secondary.weight);
    
    // 2. Convert mu to kg (CRITICAL STEP)
    double muKg = muAmu * 1.660539e-27;
    
    // 3. Estimate Force Constant (k) in N/m 
    // We scale based on Electronegativity as a proxy for bond strength
    double k = 500.0 * (1.0 + Math.abs(central.en - secondary.en) / 2.0);
    
    // 4. Calculate Frequency in Hz, then convert to THz
    // 1 THz = 10^12 Hz
    double freqHz = (1.0 / (2.0 * Math.PI)) * Math.sqrt(k / muKg);
    double freqTHz = freqHz / 1.0e12;

    System.out.printf("Vibration Frequency: %.2f THz\n", freqTHz);
}
            }

            if (totalAtoms > 1) {
                int attached = totalAtoms - counts.get(central.symbol); 
                int lp = (central.valenceElectrons - (attached)) / 2;
                if (lp < 0) lp = 0;
                int sn = attached + lp;
                System.out.println("Central Atom: " + central.name);
                System.out.println("Steric Number: " + sn + " (Atoms: " + attached + ", Lone Pairs: " + lp + ")");
                predictGeometry(sn, lp);
            }

            System.out.print("\nRun Digital Lab Simulation (MD) for this molecule? (y/n): ");
            if (sc.next().equalsIgnoreCase("y")) {
                runLabSimulation(counts, totalWeight);
            }
            
            System.out.print("\nAnalyze another? (y/n): ");
            choice = sc.next();
        }
    }

    public static void runLabSimulation(Map<String, Integer> counts, double molWeight) {
        List<Particle> particles = new ArrayList<>();
        double sigma = 3.5, epsilon = 0.2, targetTemp = 200.0;
        String molName = "Generic Compound";

        // Heuristic Intelligence Layer
        boolean hBonding = counts.containsKey("H") && (counts.containsKey("N") || counts.containsKey("O") || counts.containsKey("F"));
        if (counts.getOrDefault("C",0)==1 && counts.getOrDefault("H",0)==4) { epsilon=0.148; sigma=3.73; targetTemp=111.6; molName="Methane (CH4)"; }
        else if (counts.getOrDefault("O",0)==1 && counts.getOrDefault("H",0)==2) { epsilon=0.650; sigma=3.16; targetTemp=373.15; molName="Water (H2O)"; }
        else {
            epsilon = hBonding ? 0.45 : 0.18;
            targetTemp = (molWeight * 1.2) + (hBonding ? 180 : 60);
        }

        System.out.println("\n--- STABILIZED DIGITAL LAB: " + molName + " ---");
        int numMols = 40; double dt = 0.5;
        int side = (int) Math.ceil(Math.pow(numMols, 1.0/3.0));
        for (int i = 0; i < numMols; i++) particles.add(new Particle((i%side)*6.0, ((i/side)%side)*6.0, (i/(side*side))*6.0));

        for (int step = 0; step < 1001; step++) {
            for (Particle p : particles) { p.fx = p.fy = p.fz = 0; }
            for (int i = 0; i < particles.size(); i++) {
                for (int j = i + 1; j < particles.size(); j++) {
                    Particle p1 = particles.get(i), p2 = particles.get(j);
                    double dx = p1.x-p2.x, dy = p1.y-p2.y, dz = p1.z-p2.z;
                    double r2 = dx*dx + dy*dy + dz*dz + 0.5;
                    double r = Math.sqrt(r2);
                    if (r < 12.0) {
                        double r6 = Math.pow(sigma / r, 6);
                        double fS = (48 * epsilon / r2) * (r6 * r6 - 0.5 * r6);
                        if (fS > 80) fS = 80;
                        p1.fx += fS * dx; p1.fy += fS * dy; p1.fz += fS * dz;
                        p2.fx -= fS * dx; p2.fy -= fS * dy; p2.fz -= fS * dz;
                    }
                }
            }
            double currentKE = 0;
            for (Particle p : particles) {
                double acc = 0.0004184 / molWeight;
                p.vx += p.fx * acc * dt; p.vy += p.fy * acc * dt; p.vz += p.fz * acc * dt;
                p.vx = Math.max(-1.5, Math.min(1.5, p.vx));
                p.x += p.vx * dt; p.y += p.vy * dt; p.z += p.vz * dt;
                currentKE += (p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
            }
            double currentT = (currentKE * molWeight) / (3 * numMols * 0.008314);
            if (step % 10 == 0) {
                double scale = Math.sqrt(1.0 + 0.1 * (targetTemp / (currentT + 0.1) - 1.0));
                for (Particle p : particles) { p.vx *= scale; p.vy *= scale; p.vz *= scale; }
            }
            if (step % 200 == 0) System.out.printf("Step %d | Temp: %.2f K | Status: %s\n", step, currentT, (currentT >= targetTemp - 5 ? "STABILIZED" : "HEATING"));
        }
    }

    public static void predictPhysicalState(Map<String, Integer> counts, double weight) {
        boolean hasHBond = counts.containsKey("H") && (counts.containsKey("N") || counts.containsKey("O") || counts.containsKey("F"));
        double estimatedBP = (weight * 0.5) - 50; 
        if (hasHBond) {
            int oCount = counts.getOrDefault("O", 0);
            int nCount = counts.getOrDefault("N", 0);
            estimatedBP += (oCount * 145) + (nCount * 90); 
        }
        System.out.println("\n--- Physical State Prediction ---");
        System.out.printf("Estimated Boiling Point: ~%.1f°C\n", estimatedBP);
        if (estimatedBP < 25) System.out.println("Predicted Phase: Gas (at 25°C)");
        else if (estimatedBP < 150) System.out.println("Predicted Phase: Liquid (at 25°C)");
        else System.out.println("Predicted Phase: Solid (at 25°C)");
    }

    public static String generateFormula(Map<String, Integer> counts) {
        StringBuilder sb = new StringBuilder();
        ArrayList<String> sorted = new ArrayList<>(counts.keySet());
        Collections.sort(sorted, (a, b) -> a.equals("C") ? -1 : b.equals("C") ? 1 : a.equals("H") ? -1 : b.equals("H") ? 1 : a.compareTo(b));
        for (String s : sorted) sb.append(s).append(counts.get(s) > 1 ? counts.get(s) : "");
        return sb.toString();
    }

    public static String generateSMILES(Atom central, Map<String, Integer> counts) {
        StringBuilder sm = new StringBuilder(central.symbol);
        counts.forEach((k, v) -> {
            if (!k.equals(central.symbol)) {
                for (int i = 0; i < v; i++) if (!k.equals("H")) sm.append("(").append(k).append(")");
            }
        });
        return sm.toString();
    }

    public static void predictGeometry(int sn, int lp) {
        System.out.print("Molecular Geometry: ");
        if (sn == 2) System.out.println("Linear (180 degrees)");
        else if (sn == 3) System.out.println(lp == 0 ? "Trigonal Planar (120 degrees)" : "Bent (~118 degrees)");
        else if (sn == 4) {
            if (lp == 0) System.out.println("Tetrahedral (109.5 degrees)");
            else if (lp == 1) System.out.println("Trigonal Pyramidal (~107 degrees)");
            else System.out.println("Bent (~104.5 degrees)");
        } else System.out.println("Complex/Advanced Geometry");
    }
}
