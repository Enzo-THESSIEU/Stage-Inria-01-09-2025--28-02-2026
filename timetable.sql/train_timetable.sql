USE master;
GO
DROP TABLE dbo.Timetable;

USE darp;
GO

CREATE TABLE dbo.Timetable (
    TrainId        INT IDENTITY(1,1) NOT NULL PRIMARY KEY,
    TripName       NVARCHAR(200)      NOT NULL,
    DepartureTime  DATETIME2(0)       NOT NULL,
    ArrivalTime    DATETIME2(0)       NOT NULL,

    CONSTRAINT CK_Timetable_TimeOrder CHECK (ArrivalTime > DepartureTime)
);

-- Optional: search + ordering index
CREATE INDEX IX_Timetable_DepartureTime ON dbo.Timetable(DepartureTime);
CREATE INDEX IX_Timetable_TripName ON dbo.Timetable(TripName);
