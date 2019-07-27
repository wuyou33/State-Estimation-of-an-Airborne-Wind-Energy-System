# FehlerLog

## Alle Filter:
- wenn delta_i auf +64° dann sind die Ergebnisse besser als wenn delta_i = -64° ist. Wieso?
	- laut sprungantwort ist die dynamic bei delta_i=+64° auch besser
	- wahrscheinlich durch die nichtlinearität!
## euler ekf ahrs:

- wraptopi wurde gemacht
- varianzen wurden verringert 
- variable x = obj.EKF.x eingesetzt, anstatt direkt obj.EKF.x. Macht einen Unterschied, da in den Kalman filter obj.EKF.x verändert wird
- Sporadisch sind die Winkel im EKF komisch, obwohl nichts geändert wurde (circle = 1). Vermutung: Winkel wird am Anfang nicht richtig geschätzt (da am Anfang die Drohne nicht statisch ist, sondern sich schon bewegt) und als Initialisierung dem Kalmanfilter übergeben. Lösung: Initialisierungswerte auf einen gleichbleibenden Wert setzen z.B. gamma_i0 = [0;0;0]

  ## quaternion ekf ahrs:
- P0 richtig wählen macht sehr viel aus. Da dadurch nicht davon ausgegangen wird, dass die Anfangswerte richtig sind und der Filter somit einschwingt. Ansonsten kann es dazu kommen, dass der Filter instabil wird
- vergessen inclination in magnetometer mit einzurechnen

## eulerEKFAHRS

-  wenn man bei der eulerEKFAHRS Simulation von Markus GPS_valid=1 in Zeile 334 von perform_state_estimation__ahrs_ekf_euler.m einfügt, dann passen die Ergebnisse mit denen von meiner Simulation zusammen. (Noch nicht verstanden, warum dies so viel schlechtere Ergebnisse liefert)

## mahonyExplicitComplementaryFilter

- unterschiedliche Systemdynamik durch hinzufügen von magnetischer inklination/deklination. 
	- Ist das richtig?

## MEKF:

- Bei einem Winkel von Psi=100° ist der Filter unstabil geworden.
  - Lösung: Lag daran, dass die Quaternionmultiplikation in der falschen Reihenfolge durchgeführt wurde. (siehe commit 154b86c1747534416a45ba6e0b1d37626e108741)

## MEKFDST:

- wenn Magnetometer dauernd aus, dann ist der Filter instabil (QuaternionEKFDST hat dieses Problem nicht!)
  - Lösung:

## MEKFDST_Mag:

- wenn delta_d und delta_i größer 0 sind und man anstatt den Beschleunigungswerten die Beschleuningungswerte aus dem Geschätzten Winkel berechnet, ist der Filter instabil!
	- lag daran, dass die Anfangswerte komplett falsch waren und die Kalman Covarianzmatrix P0 mit zu geringen Werte bestückt wurde.

## euler INS

- braucht sehr lange zum einschwingen
	- Änderungen:
		- Bias Acc hinzugefügt (nichts gebracht, wieder entfernt)
		-sigma_2.a_b wurde durch TA dividiert in der initialisierung, da es jetzt ein prozessrauschen ist. (nichts gebracht)
		- rauschen ausschalten (bringt nichts)
		- inclination ausschalten (passt) --> Fehler beim rausfiltern
	- Lösung: P0 des gyro bias war viel zu klein gewählt (Bias wurde mit [1;2;3] simuliert)! --> größer machen

## euler INS DST
- braucht sehr lange zum einschwingen
	- Änderungen: 
		- nicht alle sensoren waren auf available (war nur ein bug, hat aber keine Änderungen gebracht)
		- P0 der Position auf 150 gesetzt, da nicht im Ursprung gestartet wird.
	- Lösung: entweder P0 Teil wo die Varianz der Position ist erhöhen, oder vorher die ersten Messwerte als Initialisierungen verwenden.

## MEKF INS
- Schätzt falsch
- bei sehr langen simulationen mit phi(t) =0, theta(t) = A*sin(w*t), psi(t) = 0 ist der filter instabil
	- getestet
		- incDec = 0, disable_noise=1, bias_gyro = 0: instabil
		- disable_noise = 0, incDec = 0: instabil
		- bias_gyro = 1, incDec=0, disable_noise =1: instabil
	- Lag daran, dass updateAttitude() nicht für den INS angepasst war, sondern noch für den AHRS
- Wenn der LineAngle Sensor ausgeschalten ist, von Beginn an, schwingt das system

## MEKF INSDST
- schwingt anders ein als MEKFINS
	- im correction step wird position und Geschwindigkeit unterschiedlich verändert, aber die anderen Parameter nicht
	- Lösung, in einem der measurements wurde vergessen eine diagonal matrix aus den Varianzen zu erzeugen.

## MEKF INSDST_tetherSag
- braucht sehr Lange zum Einschwingen, aber INSMEKFDST braucht das nicht
	- measurement variance war viel zu niedrig im Filter eingestellt
		- hat nichts gebracht
	- Lösung: multiplikation mit TA in der Variance wurde gelöscht, muss dort sein, da sie sonst im Filter zu groß wird, da da durch TA dividiert wird.

## quaternionEKFDifferentSampleTime
- schwingt wenn TA_magnetometer 15Hz ist, und die Frequenz des Filters bei 100Hz liegt