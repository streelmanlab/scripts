#! /usr/bin/awk -f

{OFS="\t"} 
{
	new_str1="./."
	new_str2="./."
	new_str3="./."
	new_str4="./."
	
	if (substr($10,1,3) != "./." || substr($11,1,3) != "./.") {
		new_str1=""
		if (substr($10,1,1) == "0" || substr($11,1,1) == "0") { 
			new_str1=new_str1"0" 
		} else { new_str1=new_str1"1" }
		if (substr($10,3,1) == "0" && substr($11,3,1) == "0") { 
			new_str1=new_str1"/0" 
		} else { new_str1=new_str1"/1" }

                if (substr($10,1,3) == "./.") {
                        new_str1=substr($11,1,3)
                }
                if (substr($11,1,3) == "./.") {
                        new_str1=substr($10,1,3)
                }
	}
	
        if (substr($12,1,3) != "./." || substr($13,1,3) != "./.") {
                new_str2=""
                if (substr($12,1,1) == "0" || substr($13,1,1) == "0") {
                        new_str2=new_str2"0"
                } else { new_str2=new_str2"1" }
                if (substr($12,3,1) == "0" && substr($13,3,1) == "0") {
                        new_str2=new_str2"/0"
                } else { new_str2=new_str2"/1" }

                if (substr($12,1,3) == "./.") {
                        new_str2=substr($13,1,3)
                }
                if (substr($13,1,3) == "./.") {
                        new_str2=substr($12,1,3)
                }
        }
	
        if (substr($14,1,3) != "./." || substr($15,1,3) != "./.") {
                new_str3=""
                if (substr($14,1,1) == "0" || substr($15,1,1) == "0") {
                        new_str3=new_str3"0"
                } else { new_str3=new_str3"1" }
                if (substr($14,3,1) == "0" && substr($15,3,1) == "0") {
                        new_str3=new_str3"/0"
                } else { new_str3=new_str3"/1" }
                
		if (substr($14,1,3) == "./.") {
                        new_str3=substr($15,1,3)
                }
                if (substr($15,1,3) == "./.") {
                        new_str3=substr($14,1,3)
                }
        }	
	
        if (substr($16,1,3) != "./." || substr($17,1,3) != "./.") {
                new_str4=""
                if (substr($16,1,1) == "0" || substr($17,1,1) == "0") {
                        new_str4=new_str4"0"
                } else { new_str4=new_str4"1" }
                if (substr($16,3,1) == "0" && substr($17,3,1) == "0") {
                        new_str4=new_str4"/0"
                } else { new_str4=new_str4"/1" }
                
		if (substr($16,1,3) == "./.") {
                        new_str4=substr($17,1,3)
                }
                if (substr($17,1,3) == "./.") {
                        new_str4=substr($16,1,3)
                }
        }
	
	if (new_str1 != "./." && new_str2 != "./." && new_str3 != "./." && new_str4 != "./.") {
		if (new_str1 != new_str2 || new_str1 != new_str3 || new_str1 != new_str4) {
			print $1, $2, $3, $4, $5, $6, $7, $8, $9, new_str1, new_str2, new_str3, new_str4
		}	
	}
}
