# Άσκηση 2
## Παράλληλα & Κατανεμημένα Συστήματα Υπολογιστών
## 21 Νοεμβρίου 2018

Να υλοποιήσετε σε MPI κατανεμημένο αλγόριθμο κατασκευής της δενδρικής
δομή vantage-point για ένα σύνολο δεδομένων X . Η δομή χρησιμοποιείται
στην εύρεση των k κοντινότερων γειτόνων (k nearest neighbor search –
kNN) των σημείων του vantage-point tree.

Ακολουθήστε τα παρακάτω βήματα στο πρόγραμμά σας:

1. Μετατρέψτε τον κώδικα2 ώστε να δουλεύει i) με αποστάσεις σε float
   και ii) με χωριστά υποσύνολα διαδικασιών που είναι ίσα σε μέγεθος
   και πάντα δύναμη του 2 (δείτε παρακάτω πως χρησιμοποιείται).
2. Δημιουργείστε p διεργασίες (processes) με την εκκίνηση του
   προγράμματος:
* Κάθε διεργασία έχει ένα N/p σημεία από ένα σετ δεδομένων X ∈ RN×D.
* επαναλάβετε τα παρακάτω l φορές.
* Ο αριθμός της επανάληψης προσδιορίζει σε ποια υποομάδα κάθε
διαδιακασία ανήκει, κι ποια διεργασία είναι ο “αρχηγός” για την
τρέχουσα επανάληψη. Δηλαδή, στη πρώτη επανάληψη, αρχηγός είναι η 0 για
την ομάδα 0 : p − 1, στη δεύτερη επανάληψη, η 0 είναι ο αρχηγός της
ομάδας 0 : p/2 − 1 και η p/2 της ομάδας p/2 : p − 1, στην τρίτη
επανάληψη, αρχηγοί είναι οι 0, p/4, p/2 και 3p/4 στις αντίστοιχες
υποομάδες, κ.ο.κ.
* Ο αρχηγός ανακοινώνει το vantage point που διαλέγει τυχαία από τα
  σημεία του σε όλους της ομάδας του. Η κάθε διεργασία μετράει την
  απόσταση των σημείων της από το vantage point.
* Υπολογίζεται η διάμεση απόσταση για τα στοιχεία κάθε υποομάδας με
  χρήση του κώδικα (1).
* Κάθε διεργασία χωρίζει και μετρά τα σημεία με αποστάσεις που είναι
μικρότερες ή ίσες και μεγαλύτερες από την διάμεση απόσταση και
ανακοινώνουν το πλήθος τους στον αρχηγό της υποομάδας.
* Στην συνέχεια τα στοιχεία αναδιανέμονται στην υποομάδα με την
βοήθεια του αρχηγού, ώστε οι πρώτες μισές διεργασίες να έχουν τα
στοιχεία με απόσταση μικρότερη ή ίση της διαμέσου από το vantage point
και οι υπόλοιπες με απόσταση μεγαλύτερη της διαμέσου, αντίστοιχα.
3. Υπολογίστε όλο το vantage-point tree. Μετά από l επαναλήψεις, οι
   διαδικασίες δουλεύουν μόνο με τοπικά δεδομένα χωρίς να χρειάζονται
   επικοινωνίες. Ποιο είναι το l;
4. Στο τέλος, όλες οι διεργασίες πληροφορούνται τα vantage point από
   ολόκληρο το κοινό δέντρο και φυσικά γνωρίζουν ποιος έχει ποιο
   υποδένδρο.
5. Υπολογείστε το all-k-NN.
* Υπολογίστε τους k κοντινότερους γείτονες όλων των σημείων, με k =
  2[1:8].
* Όσασημείαέχουντονkγείτοναπιομακρινόαπότηναπόστασηπουχωρίζειτουποδένδρο,ομα-
δοποιούνται και στέλνονται στον ιδιοκτήτη του υποδένδρου. Στη
συνέχεια, η διεργασία αναμέ- νει τις απαντήσεις.

